/**
 * 
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.text.ParseException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.io.LineByLineFileReader;
import edu.mit.csail.cgs.utils.io.Points2RegionsConverter;

/**
 * @author rca
 * This class parses output files from a variety of common motif discovery
 * algorithms to create WeightMatrix objects, and can also parse raw PWMs in a
 * variety of formats.
 */
public class PWMParser {
  
  private static Logger logger = Logger.getLogger(PWMParser.class);
  
  /**
   * Parse a PWM in which the rows are nucleotides and the columns are positions.
   * This method only checks that each row has the expected quantity of number
   * tokens, but doesn't make any assumptions about the format of the PWM
   * (freq, counts, log-odds, etc...). Thus each row may contain other
   * non-number tokens such as the nucleotide letter, etc...
   * @param length the number of positions in the pwm
   * @param aRow the row for A
   * @param cRow the row for C
   * @param gRow the row for G
   * @param tRow the row for T
   * @return a 4 x length DenseDoubleMatrix2D containing the parsed values
   * @throws IOException if an IO error occurs tokenizing one of the lines
   * @throws ParseException if the quantity of number tokens is not <i>length</i>
   */
  public static DenseDoubleMatrix2D parsePWM(int length, String aRow, String cRow, String gRow, String tRow) throws IOException, ParseException {
    DenseDoubleMatrix2D rawPWM = new DenseDoubleMatrix2D(4, length);
    String[] rows = new String[] {aRow, cRow, gRow, tRow};
    for (int i = 0; i < rows.length; i++) {
      String row = rows[i];
      StreamTokenizer tokenizer = new StreamTokenizer(new StringReader(row));      
      int valIndex = 0;
      boolean done = false;
      while (!done && (tokenizer.nextToken() != StreamTokenizer.TT_EOF)) {
        if (tokenizer.ttype == StreamTokenizer.TT_NUMBER) {
          rawPWM.setQuick(i, valIndex, tokenizer.nval);
          valIndex++;
          done = (valIndex == length); 
        }
      }      
      //check that no more number tokens are available
      while (tokenizer.ttype != StreamTokenizer.TT_EOF) {
        if (tokenizer.nextToken() == StreamTokenizer.TT_NUMBER) {
          valIndex++;
        }        
      }
      //check whether the right number of tokens were parsed
      if (valIndex != length) {
        throw new ParseException("Expected " + length + " numbers, but parsed " + valIndex + " from:\n" + row, row.length());
      }
    }
    
    return rawPWM;
  }
  
  
  /**
   * Parse a PWM in which the columns are nucleotides and the rows are positions.
   * This method will attempt to handle rows that include a row number as the
   * first number token. This method only checks that each column has the four 
   * number tokens (and a possible row number), but doesn't make any 
   * assumptions about the format of the PWM (freq, counts, log-odds, etc...). 
   * Thus each row may contain other non-number tokens such as the nucleotide 
   * letter, etc...
   * @param posRows An array of Strings that are the consecutive rows of the PWM
   * @return a 4 x length DenseDoubleMatrix2D containing the parsed values
   * @throws IOException if an IO error occurs tokenizing one of the lines
   * @throws ParseException if the quantity of number tokens is not 4, or 5 
   * (including a row number)
   */
  public static DenseDoubleMatrix2D parsePWM(String[] posRows) throws IOException, ParseException {
    DenseDoubleMatrix2D rawPWM = new DenseDoubleMatrix2D(4, posRows.length);
    boolean hasRowNumbers = false;
    boolean zeroBased = false;
    for (int i = 0; i < posRows.length; i++) {
      String row = posRows[i];
      StreamTokenizer tokenizer = new StreamTokenizer(new StringReader(row));      
      int valIndex = 0;
      boolean done = false;   
      boolean skippedHeader = false;
      while (!done && (tokenizer.nextToken() != StreamTokenizer.TT_EOF)) {
        if (tokenizer.ttype == StreamTokenizer.TT_NUMBER) {
          //handle row numbers
          if (hasRowNumbers && (valIndex == 0) && !skippedHeader) {
            if (zeroBased && (tokenizer.nval == i)) {
              skippedHeader = true;
            }
            else if (!zeroBased && (tokenizer.nval == i+1)) {
              skippedHeader = true;
            }
            else {
              throw new ParseException("Unknown PWM format, or incorrect row number:\n" + row, 0);
            }
          }
          else {
            rawPWM.setQuick(valIndex, i, tokenizer.nval);
            valIndex++;
            done = (valIndex == 4);
          }
        }
      }      
      //check that no more number tokens are available
      double temp = Double.NaN;
      while (tokenizer.ttype != StreamTokenizer.TT_EOF) {
        if (tokenizer.nextToken() == StreamTokenizer.TT_NUMBER) {
          temp = tokenizer.nval;
          valIndex++;
        }        
      }
      //check whether the right number of tokens were parsed
      if (valIndex != 4) {
        //If the first column may be a row number try to keep going
        if ((valIndex == 5) && (i == 0) 
            && ((rawPWM.getQuick(0, 0) == 0) || (rawPWM.getQuick(0, 0) == 1))) {
          if (rawPWM.getQuick(0, 0) == 0) {
            zeroBased = true;
          }
          else {
            zeroBased = false;
          }
          //adjust the matrix to correct for the row number
          for (int j = 0; j < 3; j++) {
            rawPWM.setQuick(j, 0, rawPWM.getQuick(j+1, 0));            
          }
          rawPWM.setQuick(3, 0, temp);
          hasRowNumbers = true;
        }
        else {
          throw new ParseException("Expected 4 numbers, but parsed " + valIndex + " from:\n" + row, row.length());
        }
      }
    }
    
    return rawPWM;
  }
  
  
  /**
   * 
   * @param priorityFile
   * @return
   * @throws IOException
   */
  public static WeightMatrix parsePriorityBestOutput(String priorityFile) throws IOException, ParseException {
    LineByLineFileReader lblfr = null;
    try {
      lblfr = new LineByLineFileReader();
      lblfr.openFile(priorityFile);
    
      char[] rowOrder = { 'A', 'C', 'G', 'T' };
      String tfString = "Transcription factor:";
      
      String tfname = null;
      DenseDoubleMatrix2D rawPWM = null;
      int pwmLength = -1;

      boolean done = false;
      String currLine = lblfr.readLine();      
      while (!done && (currLine != null)) {
        //read the tf name
        if (currLine.startsWith(tfString)) {
          tfname = currLine.substring(tfString.length()).trim();
        }
        else if (currLine.startsWith("Phi:")) {
          //read the header row
          currLine = lblfr.readLine().trim();          
          pwmLength = Integer.valueOf(currLine.substring(currLine.lastIndexOf(" ") + 1));
          
          //read the pwm
          String aRow = lblfr.readLine();
          String cRow = lblfr.readLine();
          String gRow = lblfr.readLine();
          String tRow = lblfr.readLine();
          rawPWM = PWMParser.parsePWM(pwmLength, aRow, cRow, gRow, tRow);
          done = true;
        }
        currLine = lblfr.readLine();
      }

      WeightMatrix wm = new WeightMatrix(pwmLength);
      wm.name = tfname;
      for (int i = 0; i < rowOrder.length; i++) {
        for (int j = 0; j < pwmLength; j++) {
          wm.matrix[j][rowOrder[i]] = (float)rawPWM.getQuick(i, j);
        }
      }
      return wm;      
    }
    finally {
      if (lblfr != null) {
        lblfr.closeFile();
      }
    }
  }
  
//  public static List<WeightMatrix> parsePriorityTrialOutput(String priorityFile) {
//    
//  }
//  
  
  
  /////////////////////////////////////////////////////////////////////////////
  // old code from WeightMatrixImport
  /////////////////////////////////////////////////////////////////////////////
  
  
  public static LinkedList<WeightMatrix> readTamoMatrices(String wmfile) throws IOException {
    int[] indices = { 'A', 'C', 'G', 'T' };
    LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    BufferedReader br = new BufferedReader(new FileReader(new File(wmfile)));
    String line;

    String name = null, version = null;
    WeightMatrix matrix = null;
    Vector<float[]> array = null;

    while ((line = br.readLine()) != null) {
      line = line.trim();
      if (line.length() > 0) {
        if (name == null) {
          int index = line.indexOf("\\s+");
          name = line.substring(0, index);
          version = line.substring(index, line.length()).trim();
          array = new Vector<float[]>();
        }
        else {
          String[] sarray = line.split("\\s+");
          float[] col = new float[4];
          for (int i = 0; i < col.length; i++) {
            col[i] = Float.parseFloat(sarray[i]);
          }
          array.add(col);
        }
      }
      else {
        if (name != null) {
          matrix = new WeightMatrix(array.size());
          matrix.name = name;
          matrix.version = version;
          matrix.type = "TAMO";
          for (int i = 0; i < array.size(); i++) {
            for (int j = 0; j < indices.length; j++) {
              matrix.matrix[i][indices[j]] = array.get(i)[j];
            }
          }
          matrices.add(matrix);
          // System.err.println("Added \"" + matrix.name + "\"");
        }
        array = null;
        name = null;
        version = null;
        matrix = null;
      }
    }

    if (name != null) {
      matrix = new WeightMatrix(array.size());
      matrix.name = name;
      matrix.version = version;
      for (int i = 0; i < array.size(); i++) {
        for (int j = 0; j < indices.length; j++) {
          matrix.matrix[i][indices[j]] = array.get(i)[j];
        }
      }
      matrices.add(matrix);
      // System.err.println("Added \"" + matrix.name + "\"");
    }

    br.close();
    return matrices;
  }


  /**
   * Parses weight matrices in JASPAR format. These needs the files in, eg,
   * psrg/datasets/jaspar-10_12_09/all_data/FlatFileDir. The filename should be
   * the name of the matrix_list.txt file; this method will then figure out the
   * names of the individual matrix files from that.
   * 
   * The species in matrix_list.txt are the IDs from the NCBI taxonomy database.
   * This method has a hard-coded reference to a copy of that DB in AFS since we
   * haven't loaded it anywhere else yet.
   * 
   */
  public static List<WeightMatrix> readJASPARFreqMatrices(String wmfile, String wmversion) throws IOException {
    LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    Map<Integer, String> speciesmap = new HashMap<Integer, String>();
    String taxofile = "/afs/csail.mit.edu/group/psrg/datasets/ncbi_taxonomy_nov_09/id_to_name.tsv";
    BufferedReader br = new BufferedReader(new FileReader(new File(taxofile)));
    String line = null;
    while ((line = br.readLine()) != null) {
      String pieces[] = line.split("\\t");
      speciesmap.put(Integer.parseInt(pieces[0]), pieces[1]);
    }
    br.close();

    File inputfile = new File(wmfile);
    br = new BufferedReader(new FileReader(inputfile));
    String dirname = inputfile.getParent();
    Pattern specpatt = Pattern.compile("species \"(\\d+)\"");
    Pattern grouppatt = Pattern.compile("tax_group \"(\\w+)\"");
    while ((line = br.readLine()) != null) {
      String pieces[] = line.split("\\t");
      String matrixfile = dirname + "/" + pieces[0] + ".pfm";
      WeightMatrix matrix = null;
      try {
        matrix = readJASPARFile(matrixfile);
        if (matrix == null) {
          continue;
        }
      }
      catch (IOException e) {
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
        }
        catch (NotFoundException e) {
          System.err.println("Couldn't find species " + specname + " from " + taxospecid + " in " + line);
          continue;
          // ignore it and move on
        }
      }
      else {
        matcher = grouppatt.matcher(pieces[4]);
        if (matcher.find() && matcher.group(1).equals("mammals")) {
          try {
            matrix.speciesid = (new Organism("Mus musculus")).getDBID();
          }
          catch (NotFoundException e) {
            continue; // shouldn't happen unless I typoed the name above
          }
        }
        else {
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
    while ((line = br.readLine()) != null) {
      line = line.trim();
      if (line.length() > 0) {
        String pieces[] = line.split("\\s+");
        if (matrix == null) {
          matrix = new WeightMatrix(pieces.length);
        }
        for (int i = 0; i < pieces.length; i++) {
          matrix.matrix[i][indices[lineno]] = Float.parseFloat(pieces[i]);
        }
        lineno++;
      }
    }
    br.close();
    matrix.normalizeFrequencies();
    return matrix;
  }
    public static WeightMatrix readUniProbeFile(String fname) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(new File(fname)));
        String line = null;        
        WeightMatrix matrix = null;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) {
                String pieces[] = line.split("\\s+");
                if (pieces[0].matches("[ACTG]:")) {
                    int matrixlen = pieces.length - 1;
                    if (matrix == null) {
                        matrix = new WeightMatrix(matrixlen);
                    }
                    for (int i = 1; i < pieces.length; i++) {
                        matrix.matrix[i-1][pieces[0].charAt(0)] = Float.parseFloat(pieces[i]);
                    }
                }
            }
        }
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
              version = v_string;
            }
          }
          // initialize id and accession if they're still null
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


  /**
   * Reads a matrix from a text-file. The file is assumed to have one sequence
   * per line, and all lines must be the same length. The WeightMatrix returned
   * has the frequencies of the bases at each position given.
   * 
   * @param fname
   * @return
   * @throws IOException
   * @throws ParseException
   */
  public static WeightMatrix readAlignedSequenceMatrix(String fname) throws IOException, ParseException {
    LinkedList<String> strings = new LinkedList<String>();
    String line;
    BufferedReader br = new BufferedReader(new FileReader(new File(fname)));
    while ((line = br.readLine()) != null) {
      line = line.trim();
      if (line.length() > 0) {
        strings.addLast(line);
      }
    }
    br.close();
    return WeightMatrixImport.buildAlignedSequenceMatrix(strings);
  }


  /*
   * reads a weight matrix from a file. Takes the filename as input. The file
   * must contain the log-odds matrix part of the TAMO output. Anythine else is
   * ignored. Returns a WeightMatrix filled in with whatever information is
   * available.
   */
  public static WeightMatrix readTamoMatrix(String fname) throws ParseException, FileNotFoundException, IOException {
    File wmfile = new File(fname);
    if (!wmfile.exists()) {
      throw new FileNotFoundException("Can't find file " + fname);
    }
    BufferedReader reader = new BufferedReader(new FileReader(wmfile));
    String line = reader.readLine();
    if (line == null || line.length() < 1) {
      throw new IOException("Got a null or empty line in " + fname + " when looking to skip header");
    }
    else {
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
        throw new IOException("Got a null line in " + fname + " when looking for matrix line " + i);
      }
      int index = -1;
      String[] values = line.split("\\s+");
      if (values[0].equals("#A")) {
        index = 'A';
      }
      else if (values[0].equals("#C")) {
        index = 'C';
      }
      else if (values[0].equals("#T")) {
        index = 'T';
      }
      else if (values[0].equals("#G")) {
        index = 'G';
      }
      else {
        throw new ParseException("Can't parse line " + line, 0);
      }
      for (int j = 1; j < values.length; j++) {
        results.matrix[j - 1][index] = Float.parseFloat(values[j]);
      }
    }
    return results;
  }


  /*
   * reads a weight matrix from a file. Takes the filename as input. The file
   * must contain the log-odds matrix part of the MEMO output. Anythine else is
   * ignored. Returns a WeightMatrix filled in with whatever information is
   * available.
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
    line = line.replaceFirst("^.*w=\\s*", "");
    line = line.replaceFirst("\\s*n=.*", "");
    int length = Integer.parseInt(line);
    WeightMatrix results = new WeightMatrix(length);
    for (int i = 0; i < length; i++) {
      line = reader.readLine().replaceFirst("^\\s*", "");
      lineno++;
      try {
        String[] pieces = line.split("\\s+");
        results.matrix[i]['A'] = Float.parseFloat(pieces[0]) / (float) 100.0;
        results.matrix[i]['C'] = Float.parseFloat(pieces[1]) / (float) 100.0;
        results.matrix[i]['G'] = Float.parseFloat(pieces[2]) / (float) 100.0;
        results.matrix[i]['T'] = Float.parseFloat(pieces[3]) / (float) 100.0;
      }
      catch (NumberFormatException ex) {
        System.err.println("At line " + lineno + ": " + line);
        ex.printStackTrace();
        throw ex;
      }
      catch (ArrayIndexOutOfBoundsException ex) {
        System.err.println("At line " + lineno + ": " + line);
        ex.printStackTrace();
        throw ex;
      }
    }
    return results;
  }

    public static List<WeightMatrix> readGimmeMotifsMatrices(String fname, String version, String type) throws FileNotFoundException, IOException {
        File wmfile = new File(fname);
        if (!wmfile.exists()) {
            throw new FileNotFoundException("Can't find file " + fname);
        }
        BufferedReader reader = new BufferedReader(new FileReader(wmfile));
        String line = null;
        int lineno = 1;
        List<WeightMatrix> output = new ArrayList<WeightMatrix>();
        ArrayList<String> lines = new ArrayList<String>();
        String name = null;
        while ((line = reader.readLine()) != null) {
            if (line.matches("^>.*")) {
                if (lines.size() > 0 && name != null) {
                    WeightMatrix m = parseGimmeMotifsFromLines(lines);
                    System.err.println("Parsed " +name);
                    m.name = name;
                    m.version = version;
                    m.type = type;
                    output.add(m);
                }

                name = line.replaceAll(">","");
                lines.clear();                
            } else {
                lines.add(line);
            }
        }
        if (lines.size() > 0 && name != null) {
            WeightMatrix m = parseGimmeMotifsFromLines(lines);
            m.name = name;
            m.version = version;
            m.type = type;
            output.add(m);            
        }

        reader.close();
        return output;
    }
    public static WeightMatrix parseGimmeMotifsFromLines(List<String> lines) {
        WeightMatrix m = new WeightMatrix(lines.size());
        for (int i = 0; i < lines.size(); i++) {
            String line[] = lines.get(i).split("\\s+");
            m.matrix[i]['A'] = Float.parseFloat(line[0]);
            m.matrix[i]['C'] = Float.parseFloat(line[1]);
            m.matrix[i]['G'] = Float.parseFloat(line[2]);
            m.matrix[i]['T'] = Float.parseFloat(line[3]);
        }
        return m;
    }

  public static void main(String[] args) {
//    int length = 6;
//    String aLine = " a  1.0000 0.9487 0.5641 0.0000 0.0000 0.7949"; 
//    String cLine = " c  0.0000 0.0513 0.2051 0.3590 1.0000 0.0000"; 
//    String gLine = " g  0.0000 0.0000 0.1538 0.2436 0.0000 0.0641";
//    String tLine = " t  0.0000 0.0000 0.0769 0.3974 0.0000 0.1410"; 

    String[] rows = new String[] {
        "01    17    22    17    43",
        "02    12    41    35    12",
        "03    34     9    56     0",
        "04    11    48    33     9",
        "05     6    80     8     6",
        "06     4    39    57     0",
        "07    91     0     3     6",
        "08     0     0     0   100",
        "09     0     0     0   100",
        "10    85     0     6     8",
        "11    45    13    16    26",
        "12     7    37    49     6",
        "13    18    24    35    24",
        "14    33    33    33     0",
        "15    25    25    25    25",
        "16    17    17    22    43"
    };
    
    try {
//      PWMParser.parsePWM(length, aLine, cLine, gLine, tLine);
      DoubleMatrix2D rawPWM = PWMParser.parsePWM(rows);
      System.out.println(rawPWM.toString());
    }
    catch (IOException ioex) {
      logger.fatal(ioex);
    }
    catch (ParseException pex) {
      logger.fatal(pex); 
    }
  }
}
