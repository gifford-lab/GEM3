package edu.mit.csail.cgs.datasets.motifs;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.sql.*;
import java.text.ParseException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.species.Organism;

/** Imports a weight matrix from a TAMO formatted file.
 *  Only reads the first WM in the file. 
 *
 * Usage:
 * java edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport --species "Saccharomyces cerevisiae" --wmname HSF1 --wmversion MacIsaac06 --wmtype TAMO --wmfile v1.HSF1.wm
*/

public class WeightMatrixImport {

  private static int MAX_MOTIF_LEN = 200;
  
    public static void main(String args[]) {
        String species = null;
        String wmname = null, wmversion = null, wmtype = null;
        String wmfile = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) { 
                species = args[++i];
                if (species.indexOf(';') != -1) {
                    String[] pieces = species.split(";");
                    species = pieces[0];
                }
            }
            if (args[i].equals("--wmname")) {
                wmname = args[++i];
                if (wmname.indexOf(';') != -1) {
                    String[] pieces = wmname.split(";");
                    wmname = pieces[0];
                    wmversion = pieces[1];
                    if (pieces.length >= 3) {
                        wmtype = pieces[2];
                    }
                }
            }
            if (args[i].equals("--wmversion")) {
                wmversion = args[++i];
            }
            if (args[i].equals("--wmtype")) {
                wmtype = args[++i];
            }
            if (args[i].equals("--") ||
                args[i].equals("--wmfile")) {
                wmfile = args[++i];
            }
        }

        //        if (species == null) {
        //            System.err.println("Must supply a --species"); System.exit(1);
        //        }
        if (wmfile == null) {
            System.err.println("Must supply a --wmfile"); System.exit(1);
        } 
        try {
            if(wmname==null) { 
                insertMultiWMFromFile(species,wmtype,wmfile, wmversion);
            } else { 
                if (wmversion == null) {
                    System.err.println("Must supply a --wmversion"); System.exit(1);
                }
                insertWMFromFile(species,wmname,wmversion,wmtype,wmfile);
            }
        } catch (SQLException ex) {
            ex.printStackTrace();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            System.err.println("Must supply a valid species and genome");
        } catch (UnknownRoleException ex ){
            ex.printStackTrace();
            System.err.println("Couldn't connect to role annotations");
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
            System.err.println("Couldn't find the input file");
        } catch (ParseException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    public static int insertMatrixIntoDB(WeightMatrix matrix) 
        throws SQLException, NotFoundException {
        
        java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
        int wmid = -1;
        PreparedStatement exists = cxn.prepareStatement("select id from weightmatrix where species = ? and name = ? and version = ? and type = ?");
        exists.setInt(1,matrix.speciesid);
        exists.setString(2,matrix.name);
        exists.setString(3,matrix.version);
        exists.setString(4,matrix.type);
        ResultSet rs = exists.executeQuery();
        if (!rs.next()) {
            rs.close();
            PreparedStatement insertwm = cxn.prepareStatement("insert into weightmatrix(id,species,name,version,type) values (weightmatrix_id.nextval,?,?,?,?)");
            insertwm.setInt(1,matrix.speciesid);
            insertwm.setString(2,matrix.name);
            insertwm.setString(3,matrix.version);
            insertwm.setString(4,matrix.type);
            //            System.err.println(String.format("Inserting %s %s %s", matrix.name, matrix.version, matrix.type));
            insertwm.execute();
            insertwm.close();
            PreparedStatement getId = cxn.prepareStatement("select weightmatrix_id.currval from dual");
            rs = getId.executeQuery();
            if (rs.next()) {
                wmid = rs.getInt(1);                
            } else {
                System.err.println("No weight matrix id");
                System.exit(1);
            }            
            rs.close();
            getId.close();
        } else {
            wmid = rs.getInt(1);
            PreparedStatement deleteold = cxn.prepareStatement("delete from weightmatrixcols where weightmatrix = ?");
            deleteold.setInt(1,wmid);
            deleteold.execute();
            deleteold.close();
        }
        rs.close();
        exists.close();
           
            
        PreparedStatement insertcol = cxn.prepareStatement("insert into weightmatrixcols(weightmatrix,position,letter,weight) values(?,?,?,?)");
        insertcol.setInt(1,wmid);
        for (int col = 0; col < matrix.length(); col++) {
//             System.err.println(String.format("  Column %d %f %f %f %f", col,
//                                              matrix.matrix[col]['A'],
//                                              matrix.matrix[col]['C'],
//                                              matrix.matrix[col]['T'],
//                                              matrix.matrix[col]['G']));
            insertcol.setInt(2,col);

            insertcol.setString(3,"A");
            insertcol.setFloat(4,matrix.matrix[col]['A']);
            insertcol.execute();
                
            insertcol.setString(3,"C");
            insertcol.setFloat(4,matrix.matrix[col]['C']);
            insertcol.execute();

            insertcol.setString(3,"T");
            insertcol.setFloat(4,matrix.matrix[col]['T']);
            insertcol.execute();

            insertcol.setString(3,"G");
            insertcol.setFloat(4,matrix.matrix[col]['G']);
            insertcol.execute();                
        }           
        insertcol.close();
        DatabaseFactory.freeConnection(cxn);
        return wmid;
    }
    
    public static LinkedList<WeightMatrix> readTamoMatrices(String wmfile) throws IOException {
        int[] indices = { 'A', 'C', 'G', 'T' };
        LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
        BufferedReader br = new BufferedReader(new FileReader(new File(wmfile)));
        String line;
        
        String name = null, version=null;
        WeightMatrix matrix = null;
        Vector<float[]> array = null;
        
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) {
                if(name == null) { 
                    int index = line.indexOf("\\s+");
                    name = line.substring(0, index);
                    version = line.substring(index, line.length()).trim();
                    array = new Vector<float[]>();
                } else { 
                    String[] sarray = line.split("\\s+");
                    float[] col = new float[4];
                    for(int i = 0; i < col.length; i++) { 
                        col[i] = Float.parseFloat(sarray[i]);
                    }
                    array.add(col);
                }
            } else { 
                if(name != null) { 
                    matrix = new WeightMatrix(array.size());
                    matrix.name=name;
                    matrix.version=version;
                    matrix.type="TAMO";
                    for(int i = 0; i < array.size(); i++) { 
                        for(int j = 0; j < indices.length; j++) {
                            matrix.matrix[i][indices[j]] = array.get(i)[j];
                        }
                    }
                    matrices.add(matrix);
                    //                    System.err.println("Added \"" + matrix.name + "\"");
                }
                array = null;
                name = null;
                version=null;
                matrix = null;
            }
        }

        if(name != null) { 
            matrix = new WeightMatrix(array.size());
            matrix.name=name;
            matrix.version=version;
            for(int i = 0; i < array.size(); i++) { 
                for(int j = 0; j < indices.length; j++) {
                    matrix.matrix[i][indices[j]] = array.get(i)[j];
                }
            }
            matrices.add(matrix);
            //            System.err.println("Added \"" + matrix.name + "\"");
        }

        br.close();
        return matrices;
    }
    /** Parses weight matrices in JASPAR format.  These needs the files in, eg, psrg/datasets/jaspar-10_12_09/all_data/FlatFileDir.  
     * The filename should be the name of the matrix_list.txt file; this method will then figure out the names
     * of the individual matrix files from that.  
     *
     * The species in matrix_list.txt are the IDs from the NCBI taxonomy database.  This method has a hard-coded
     * reference to a copy of that DB in AFS since we haven't loaded it anywhere else yet.
     *
     */
    public static List<WeightMatrix> readJASPARFreqMatrices(String wmfile, String wmversion) throws IOException {
        LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
        Map<Integer,String> speciesmap = new HashMap<Integer,String>();
        String taxofile = "/afs/csail.mit.edu/group/psrg/datasets/ncbi_taxonomy_nov_09/id_to_name.tsv";
        BufferedReader br = new BufferedReader(new FileReader(new File(taxofile)));
        String line = null;        
        while((line = br.readLine()) != null) { 
            String pieces[] = line.split("\\t");
            speciesmap.put(Integer.parseInt(pieces[0]), pieces[1]);
        }
        br.close();

        File inputfile = new File(wmfile);
        br = new BufferedReader(new FileReader(inputfile));
        String dirname = inputfile.getParent();
        Pattern specpatt = Pattern.compile("species \"(\\d+)\"");
        Pattern grouppatt = Pattern.compile("tax_group \"(\\w+)\"");
        while((line = br.readLine()) != null) { 
            String pieces[] = line.split("\\t");
            String matrixfile = dirname + "/" + pieces[0] + ".pfm";
            WeightMatrix matrix = null;
            try {
                matrix = readJASPARFile(matrixfile);
                if (matrix == null) {continue;}
            } catch (IOException e) {
                System.err.println("Couldn't read matrixfile " + e.toString());
                continue;
            }
            matrix.name = pieces[2];
            matrix.type = "JASPAR";
            matrix.version = wmversion + " " + pieces[0];
            Matcher matcher = specpatt.matcher(pieces[4]);
            if (matcher.find()) {
                int taxospecid = Integer.parseInt(matcher.group(1));
                String specname = speciesmap.get(taxospecid);
                try {
                    matrix.speciesid = (new Organism(specname)).getDBID();
                    matrix.species = specname;
                } catch (NotFoundException e) {
                    System.err.println("Couldn't find species " + specname + " from " + taxospecid + " in " + line);
                    continue;
                    // ignore it and move on
                }
            } else {
                matcher = grouppatt.matcher(pieces[4]);
                if (matcher.find() && matcher.group(1).equals("mammals")) {
                    try {
                        matrix.speciesid = (new Organism("Mus musculus")).getDBID();
                    } catch (NotFoundException e) {
                        continue;  // shouldn't happen unless I typoed the name above
                    }
                } else {
                    System.err.println("No species in " + line);
                    continue;
                }
            }
            matrices.add(matrix);
        }
        br.close();
        return matrices;
    }
    public static WeightMatrix readJASPARFile(String fname) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(new File(fname)));
        String line = null;        
        WeightMatrix matrix = null;
        int[] indices = { 'A', 'C', 'G', 'T' };
        int lineno = 0;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) {
                String pieces[] = line.split("\\s+");
                if (matrix == null) {
                    matrix = new WeightMatrix(pieces.length);                    
                }
                for (int i = 0; i< pieces.length; i++) {
                    matrix.matrix[i][indices[lineno]] = Float.parseFloat(pieces[i]);
                }
                lineno++;
            }
        }
        br.close();
        matrix.normalizeFrequencies();
        return matrix;
    }
    
    
  /**
   * Parses in weight matrices in transfac format with counts or frequencies
   * (not log odds)
   * 
   * @param wmfile
   * @param version
   * @return
   * @throws IOException
   */
  public static List<WeightMatrix> readTRANSFACFreqMatrices(String wmfile, String version) throws IOException {
    int[] indices = { 'A', 'C', 'G', 'T' };
    LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    BufferedReader br = new BufferedReader(new FileReader(new File(wmfile)));
    String line;
    WeightMatrix matrix = null;
    int motifCount = 0;
    Vector<float[]> arrays = new Vector<float[]>();

    // Read in Transfac format first
    Organism currentSpecies = null;
    String name = null, id = null, accession = null;
    Pattern speciesPattern = Pattern.compile(".*Species:.*, (.*)\\.");
    while ((line = br.readLine()) != null) {
      line = line.trim();
      if (line.length() > 0) {
        String[] pieces = line.split("\\s+");
        if (pieces[0].equals("AC")) {
          accession = pieces[1];
          // System.err.println("read AC " + accession);
        }
        else if (pieces[0].equals("NA")) {
          name = pieces[1];
          // System.err.println("read NA " + name);
        }
        else if (pieces[0].equals("ID")) {
          id = pieces[1];
          // System.err.println("read ID " + id);
          arrays.clear();
          name = null;
          accession = null;
          currentSpecies = null;
        }
        else if (pieces[0].equals("BF") && currentSpecies == null) {
          Matcher matcher = speciesPattern.matcher(line);
          if (matcher.matches()) {
            String specname = matcher.group(1);
            try {
              currentSpecies = new Organism(specname);
              // System.err.println("Got species " + specname);
            }
            catch (NotFoundException e) {
              System.err.println("Couldn't find species " + specname);
              // ignore it and move on
            }
          }
        }
        else if (pieces[0].equals("DE")) {
          name = pieces[1];
          if (pieces.length >= 3) {
            String v_string = pieces[2];
            if (pieces.length >= 4) {
              for (int v = 3; v < pieces.length; v++) {
                v_string = v_string + "," + pieces[v];
              }
            }

            if (version != null) {
              version = v_string + "," + version;
            }
            else {
              version  = v_string;
            }
          }        
          //initialize id and accession if they're still null
          if (id == null) {
            id = "";           
          }
          if (accession == null) {
            accession = "";
          }
        }
        else if (pieces[0].equals("XX")) {
          if (name != null && accession != null && id != null && arrays.size() > 0) {
            matrix = new WeightMatrix(arrays.size());
            for (int i = 0; i < arrays.size(); i++) {
              matrix.matrix[i]['A'] = arrays.get(i)[0];
              matrix.matrix[i]['C'] = arrays.get(i)[1];
              matrix.matrix[i]['G'] = arrays.get(i)[2];
              matrix.matrix[i]['T'] = arrays.get(i)[3];
            }
            matrix.name = name;
            matrix.version = version;
            if (id.length() > 0) {
              matrix.version = matrix.version + " " + id;
            }
            if (accession.length() > 0) {
              matrix.version = matrix.version + " " + accession;
            }
            
            // System.err.println("read version " + matrix.version);
            if (currentSpecies != null) {
              matrix.speciesid = currentSpecies.getDBID();
              matrix.species = currentSpecies.getName();
            }
            matrix.type = "TRANSFAC";
            matrix.normalizeFrequencies();
            matrices.add(matrix);

            // clean up to prepare to parse next pwm
            arrays.clear();
            name = null;
            id = null;
            accession = null;
            currentSpecies = null;
          }
          else {
            // System.err.println(String.format("name %s id %s species %s",name,id,currentSpecies
            // == null ? "null" : currentSpecies.toString()));
          }
        }
        else if (name != null && (pieces.length == 5 || pieces.length == 6) && Character.isDigit(pieces[0].charAt(0))) {
          // Load the matrix
          float[] fa = new float[4];
          // System.err.println("  adding matrix line");
          for (int i = 1; i <= 4; i++) {
            fa[i - 1] = Float.parseFloat(pieces[i]);
          }
          arrays.add(fa);
        }
      }
    }
    
    br.close();
    return matrices;
  }


    public static Set<Integer> insertMultiWMFromFile(String species, String wmtype, String wmfile, String wmversion) 
        throws IOException, SQLException, NotFoundException { 

        HashSet<Integer> ids = new HashSet<Integer>();
        List<WeightMatrix> matrices = null;
        if(wmtype.matches(".*TAMO.*")) { 
            int speciesid = (new Organism(species)).getDBID();
            matrices = readTamoMatrices(wmfile);
            for(WeightMatrix matrix : matrices) { 
                matrix.speciesid = speciesid;
            }
        } else if (wmtype.matches(".*TRANSFAC.*")) {
            matrices = readTRANSFACFreqMatrices(wmfile, wmversion);
        } else if (wmtype.matches(".*JASPAR.*")) {
            matrices = readJASPARFreqMatrices(wmfile, wmversion);
        } else {
            throw new NotFoundException("Unknown weight matrix type " + wmtype);
        }
            for(WeightMatrix matrix : matrices) { 
                ids.add(insertMatrixIntoDB(matrix));                
            }
        return ids;
    }

    /* reads a weight matrix from the specified file,
       inserts it into the database with the specified name and
       version, and returns its dbid.  This means that the file may contain
       only a single weight matrix*/
    public static int insertWMFromFile(String species,
                                       String wmname,
                                       String wmversion,
                                       String wmtype,
                                       String wmfile) throws SQLException, NotFoundException, UnknownRoleException, FileNotFoundException, ParseException, IOException {
        WeightMatrix matrix;
        if (wmtype.matches(".*TAMO.*")) {
            matrix = readTamoMatrix(wmfile);
        } else if (wmtype.matches(".*MEME.*")) {
            matrix = readMemeMatrix(wmfile);
        } else if (wmtype.matches(".*SEQ.*")) { 
            matrix = readAlignedSequenceMatrix(wmfile);
        } else if (wmtype.matches(".*TRANSFAC.*")) {
          //TODO add a method to read a single transfac matrix
          matrix = readTRANSFACFreqMatrices(wmfile, wmversion).get(0);
        }
        else {
            System.err.println("Didn't see a program I recognize in the type.  defaulting to reading TAMO format");
            matrix = readTamoMatrix(wmfile);
        }
        matrix.name = wmname;
        matrix.version = wmversion;
        matrix.type = wmtype;
        matrix.speciesid = (new Organism(species)).getDBID();
        return insertMatrixIntoDB(matrix);
    }
    
    /**
     * Constructs a matrix from a set of strings.  The strings must all have the same length.
     * The WeightMatrix returned has the frequencies of the bases at each position given.
     * 
     * @param strings
     * @return
     * @throws IOException
     * @throws ParseException
     */
    public static WeightMatrix buildAlignedSequenceMatrix(Collection<String> strings) throws ParseException {
        WeightMatrix wm = null;
        
        int[] counts = null;
        for(String line : strings) { 
            line = line.trim().toUpperCase();
            if(line.length() > 0) { 
                if(wm == null) { 
                    wm = new WeightMatrix(line.length());
                    counts = new int[line.length()];
                    for(int i = 0; i < wm.length(); i++) {
                        counts[i] = 0;
                        for(int j = 0; j < wm.matrix[i].length; j++) { 
                            wm.matrix[i][j] = (float)0.0;
                        }
                    }
                }
                
                if(line.length() != wm.length()) { 
                    throw new ParseException("Line \"" + line + "\" was of uneven length (" + 
                            wm.length() + ")", 0);
                }
                
                for(int i = 0; i < line.length(); i++) {
                    char c = line.charAt(i);
                    if(c != 'N' && c != '-') { 
                        wm.matrix[i][c] += (float)1.0;
                        counts[i] += 1;
                    }
                }
            } 
        }
        
        for(int i = 0; wm != null && i < wm.length(); i++) { 
            if(counts[i] > 0) { 
                for(int j = 0; j < wm.matrix[i].length; j++) {
                    wm.matrix[i][j] /= (float)counts[i];
                }
            } else { 
                for(int j = 0; j < wm.matrix[i].length; j++) { 
                    wm.matrix[i][j] = (float)1.0 / (float)(wm.matrix[i].length);
                }
            }
        }
        
        return wm;
    }

    /**
     * Reads a matrix from a text-file.  The file is assumed to have one sequence per line, 
     * and all lines must be the same length.  The WeightMatrix returned has the frequencies 
     * of the bases at each position given.
     * 
     * @param fname
     * @return
     * @throws IOException
     * @throws ParseException
     */
    public static WeightMatrix readAlignedSequenceMatrix(String fname) throws IOException,ParseException {
        LinkedList<String> strings = new LinkedList<String>();
        String line;
        BufferedReader br = new BufferedReader(new FileReader(new File(fname)));
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                strings.addLast(line);
            }
        }
        br.close();
        return buildAlignedSequenceMatrix(strings);
    }

    /* reads a weight matrix from a file.  Takes the filename as input.
       The file must contain the log-odds matrix part of the TAMO output.  Anythine
       else is ignored.
       Returns a WeightMatrix filled in with whatever information is available.
    */
    public static WeightMatrix readTamoMatrix(String fname) throws ParseException, FileNotFoundException, IOException {
        File wmfile = new File(fname);
        if (!wmfile.exists()) {
            throw new FileNotFoundException("Can't find file " + fname);
        }
        BufferedReader reader = new BufferedReader(new FileReader(wmfile));
        String line = reader.readLine();
        if (line == null ||
            line.length() < 1) {
            throw new IOException("Got a null or empty line in " + fname + " when looking to skip header");
        } else {
            System.err.println("READ " + line);
        }
        while (line.charAt(0) != '#') {
            line = reader.readLine();
            System.err.println("READ " + line);
        }
        String[] positions = line.split("\\s+");
        Integer length = Integer.parseInt(positions[positions.length - 1]) + 1;
        System.err.println("Length is " + length);
        WeightMatrix results = new WeightMatrix(length);
        for (int i = 0; i <= 3; i++) {
            line = reader.readLine();
            if (line == null) {
                throw new IOException("Got a null line in " +fname+" when looking for matrix line " + i);
            }
            int index = -1;
            String[] values = line.split("\\s+");
            if (values[0].equals("#A")) {
                index = 'A';
            } else if (values[0].equals("#C")) {
                index = 'C';
            }
            else if (values[0].equals("#T")) {
                index = 'T';
            }
            else if (values[0].equals("#G")) {
                index = 'G';
            } else {
                throw new ParseException("Can't parse line " + line,0);
            }
            for (int j = 1; j < values.length; j++) {
                results.matrix[j-1][index] = Float.parseFloat(values[j]);
            }
        }
        return results;
    }


    /* reads a weight matrix from a file.  Takes the filename as input.
       The file must contain the log-odds matrix part of the MEMO output.  Anythine
       else is ignored.
       Returns a WeightMatrix filled in with whatever information is available.
    */
    public static WeightMatrix readMemeMatrix(String fname) throws ParseException, FileNotFoundException, IOException {
        File wmfile = new File(fname);
        if (!wmfile.exists()) {
            throw new FileNotFoundException("Can't find file " + fname);
        }
        BufferedReader reader = new BufferedReader(new FileReader(wmfile));
        String line = reader.readLine();
        int lineno = 1;
        while (!line.matches(".*log-odds matrix.*")) {
            line = reader.readLine();
            lineno++;
        }
        line = line.replaceFirst("^.*w=\\s*","");
        line = line.replaceFirst("\\s*n=.*","");
        int length = Integer.parseInt(line);
        WeightMatrix results = new WeightMatrix(length);
        for (int i = 0; i < length; i++) {
            line = reader.readLine().replaceFirst("^\\s*","");
            lineno++;
            try {
                String[] pieces = line.split("\\s+");
                results.matrix[i]['A'] = Float.parseFloat(pieces[0]) / (float)100.0;
                results.matrix[i]['C'] = Float.parseFloat(pieces[1]) / (float)100.0;
                results.matrix[i]['G'] = Float.parseFloat(pieces[2]) / (float)100.0;
                results.matrix[i]['T'] = Float.parseFloat(pieces[3]) / (float)100.0;
            } catch (NumberFormatException ex) {
                System.err.println("At line " + lineno + ": "+ line);
                ex.printStackTrace();
                throw ex;
            } catch (ArrayIndexOutOfBoundsException ex) {
                System.err.println("At line " + lineno + ": "+ line);
                ex.printStackTrace();
                throw ex;
            }
        }
        return results;
    }

}

