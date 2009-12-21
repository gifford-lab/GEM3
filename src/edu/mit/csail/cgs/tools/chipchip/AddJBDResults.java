package edu.mit.csail.cgs.tools.chipchip;

import java.io.*;
import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.parsing.textfiles.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;

/**
 * Stores JBD results to the database.  gse-db-scripts/run_jbd.pl will generate a suitable command line
 * for this program (or the perl version).
 *
 * java edu.mit.csail.cgs.tools.chipchip.AddJBDResults --analysis 'jbd analysis name;analysis version' --species "$MM;mm8" \
 *    --paramfile params.txt \
 *    [--pattern chr10] \
 *    [--chromindex 3] \
 *    [--dir ./ ] \
 *    [--expt 'expt name;version;replicate'] \
 *    -- *.jbd
 *
 *   --paramfile is a set of key=value pairs that give the parameters used for the analysis.  These
 *     are stored in the db.  This file is produced by ep_cmdline.pl when the JBD command lines are generated
 *   
 *   --chromindex : file names are split on . and chromindex tells AddJBDResults which 0-indexed
 *     piece is the chromosome name.  Default value is 1
 *
 *   --pattern : regular expression.  Only jbd results files matching this pattern are loaded.
 *
 *   --dir : directory in which to look for input files.  All files with the extension '.jbd' (and matching
 *     the optional --pattern) are read
 *  
 *   --expt : the name;version;replicate of an experiment on which this analysis was based.  Can be specified
 *     more than once.  The database records this along with the values from the paramfile so that
 *     an analysis could be recreated (or rerun) based on the information in the db.
 */

public class AddJBDResults {

    private String analysisname, analysisversion, speciesname, genomename, paramfile, pattern;
    private Set<String> directories;
    private int chromindex;
    private Set<Experiment> experiments;
    private ChipChipMetadataLoader loader;
    private Organism organism;
    private Genome genome;
    private List<File> files;
    private PreparedStatement addObs;
    private java.sql.Connection cxn;

    public static void main(String args[]) {
        AddJBDResults adder = null;
        try {
            adder = new AddJBDResults();
            adder.parseArgs(args);
            adder.fillArgs();
            adder.add();
        } catch (Exception e) {
            System.err.println("Couldn't add files");
            e.printStackTrace();
            try {
                adder.cxn.rollback();
            } catch (SQLException sqlex) {
                System.err.println("Couldn't rollback");
                sqlex.printStackTrace();
            }
        }
    }

    public AddJBDResults () throws SQLException {
        loader = new ChipChipMetadataLoader();
        analysisname = null;
        analysisversion = null;
        speciesname = null;
        genomename = null;
        paramfile = null;
        pattern = null;
        directories = new HashSet<String>();
        chromindex = 1;
        experiments = new HashSet<Experiment>();
        files = new ArrayList<File>();
        cxn = edu.mit.csail.cgs.utils.database.DatabaseFactory.getConnection("chipchip");
        cxn.setAutoCommit(false);
        addObs = cxn.prepareStatement("insert into bayesresults(analysis,chromosome,position,posterior,posteriorstd,strength,strengthstd)"+
                                      " values (?,?,?,?,?,?,?)");
    }

    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--analysis")) {
                String pieces[] = args[++i].split(";");
                analysisname = pieces[0];
                analysisversion = pieces[1];
            }
            if (args[i].equals("--analysisname")) {
                analysisname = args[++i];
            }
            if (args[i].equals("--analysisversion")) {
                analysisversion = args[++i];
            }
            if (args[i].equals("--species")) {
                String pieces[] = args[++i].split(";");
                speciesname = pieces[0];
                if (pieces.length >= 1) {
                    genomename = pieces[1];
                }
            }
            if (args[i].equals("--genome")) {
                genomename = args[++i];
            }
            if (args[i].equals("--paramfile")) {
                paramfile = args[++i];
            }
            if (args[i].equals("--pattern")) {
                pattern = args[++i];
            }
            if (args[i].equals("--chromindex")) {
                chromindex = Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--dir")) {
                directories.add(args[++i]);
            }
            if (args[i].equals("--expt")) {
                String pieces[] = args[++i].split(";");
                if (pieces.length == 2) {
                    experiments.addAll(loader.loadExperiment(pieces[0],pieces[1]));
                } else if (pieces.length == 3) {
                    experiments.add(loader.loadExperiment(pieces[0],pieces[1],pieces[2]));
                } else {
                    throw new NotFoundException("Couldn't figure out experiment " + pieces);
                }
            }
            if (args[i].equals("--")) {
                for (i++; i < args.length; i++) {
                    files.add(new File(args[i]));
                }
            }
        }        
        organism = new Organism(speciesname);
        genome = organism.getGenome(genomename);        
    }
    
    /* make sure that all required args have values.  If they don't, try to guess them.
       throws InvalidArgumentException if it no suitable value for required 
       arguments can be found.  */       
    public void fillArgs() throws IllegalArgumentException {
        if (paramfile == null) {
            for (String d : directories) {
                File f= new File(d + File.separator + "bayes.params");
                if (f.exists()) {
                    paramfile = f.getPath();
                }
            }
        }
        if (paramfile == null) {
            throw new IllegalArgumentException("No value for parameters file");
        }
        File params = new File(paramfile);
        if (!params.exists()) {
            throw new IllegalArgumentException("Parameters file " + paramfile + " doesn't exist");
        }              
    }

    Map<String,String> readParamFile(String fname) throws IOException {
        HashMap<String,String> output = new HashMap<String,String>();
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line;
        while ((line = reader.readLine()) != null) {
            String[] pieces = line.split("\\t");
            output.put(pieces[0],pieces[1]);
        }
        reader.close();
        return output;
    }

    /* look for .jbd files in the directories specified in directories */
    public void expandDirectories() throws IOException {
        for (String dirname : directories) {
            File directory = new File(dirname);
            if (!directory.exists() && directory.isDirectory()) {
                continue;
            }
            expandDirectory(directory);
        }
    }
    /* recursively add files in the specified directory to
       files
    */
    private void expandDirectory(File directory) throws IOException {
        File[] dirfiles = directory.listFiles();        
        if (!directory.exists() ||
            !directory.isDirectory()) {
            return;
        }
        for (File file : dirfiles) {
            if (file.isDirectory() && !file.getName().matches("^\\..*")) {
                expandDirectory(file);
            } else {
                files.add(file);
            }
        }
    }

    /* filters files to weed out
       - files that don't exist
       - files whose name doesn't end in .jbd
       - files whose name doesn't match the pattern, if provided */
    public void filterFiles() throws IOException {
        ArrayList<File> out = new ArrayList<File>();
        for (File file : files) {
            if (!file.exists() && file.isFile()) {
                continue;
            }
            if (!file.getName().matches(".*\\.jbd$")) {
                continue;
            }
            if (pattern == null) {
                out.add(file);
                continue;
            }
            if (file.getName().matches(pattern)) {
                out.add(file);
            }
        }
        files = out;
    }

    public void addFile(File file) throws SQLException, IOException, NotFoundException {
        String[] pieces = file.getName().split("\\.");
        String chromnum = pieces[chromindex];
        int chromid = genome.getChromID(chromnum);
        addObs.setInt(2,chromid);
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        while ((line = reader.readLine()) != null) {
            pieces = line.split("\\t");
            if (pieces.length != 5) {
                throw new IOException("Too few fields in line " + line);
            }
            addObs.setInt(3,Integer.parseInt(pieces[0]));
            addObs.setDouble(4,Double.parseDouble(pieces[1]));
            addObs.setDouble(5,Double.parseDouble(pieces[2]));
            addObs.setDouble(6,Double.parseDouble(pieces[3]));
            addObs.setDouble(7,Double.parseDouble(pieces[4]));
            addObs.execute();
        }        
        reader.close();
    }

    public void add() throws SQLException, IOException, NotFoundException {
        JBDAnalysis analysis = loader.loadJBDAnalysis(organism.getDBID(),
                                                      analysisname,
                                                      analysisversion);
        Map<String,String> params = readParamFile(paramfile);        
        loader.addJBDParam(analysis,params);
        for (Experiment e : experiments) {
            loader.addJBDInput(analysis,e);
        }
        loader.mapJBDAnalysisToGenome(analysis.getDBID(), genome.getDBID());
        expandDirectories();
        filterFiles();
        addObs.setInt(1,analysis.getDBID());
        for (File file : files) {
            addFile(file);
        }
        cxn.commit();
    }

}
