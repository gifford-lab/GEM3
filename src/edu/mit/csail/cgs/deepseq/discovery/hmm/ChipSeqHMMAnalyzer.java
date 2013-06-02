/**
 * 
 */
package edu.mit.csail.cgs.deepseq.discovery.hmm;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixImport;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.utilities.ChipSeqAnalysisUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.io.parsing.PWMParser;

/**
 * @author rca
 * 
 */
public class ChipSeqHMMAnalyzer {

  private static Logger logger = Logger.getLogger(ChipSeqHMMAnalyzer.class);

  private static final String FRAG_DIST_FILENAME_KEY = "frag_file";
  private static final String MODEL_MIN_KEY = "model_min";
  private static final String MODEL_MAX_KEY = "model_max";
  private static final String SELECTED_REGION_FILENAME_KEY = "region_file";
  private static final String SELECTED_REGION_FILE_FORMAT_KEY = "region_file_format";
  private static final String MINIMUM_READS_FOR_SELECTED_REGION_KEY = "min_reads";
  private static final String MOTIF_PWM_NAME = "motif_PWM_name";
  private static final String MOTIF_PWM_VERSION = "motif_PWM_ver";
  private static final String BG_PWM_NAME = "bg_PWM_name";
  private static final String BG_PWM_VERSION = "bg_PWM_ver";
  private static final String MOTIF_PWM_FILENAME = "motif_PWM_filename";
  private static final String BG_PWM_FILENAME = "bg_PWM_filename";
  
  private static final String FRAG_DIST_FILENAME_DFLT = "oct4.shear.txt";
  private static final int MODEL_MIN_DFLT = -200;
  private static final int MODEL_MAX_DFLT = 200;
  private static final String SELECTED_REGION_FILENAME_DFLT = "selected_regions.txt";
  private static final String SELECTED_REGION_FILE_FORMAT_DFLT = "Regions";
  private static final double MINIMUM_READS_FOR_SELECTED_REGION_DFLT = 3.0;

  
  private String[] args;

  private ArrayList<Pair<DeepSeqExpt, DeepSeqExpt>> experiments = new ArrayList<Pair<DeepSeqExpt, DeepSeqExpt>>();

  private Genome genome;
  
  private int readLength=32;

  private BindingMotifHMM chipSeqModel = null;

  /**
   * This constructor will get most parameters from property file
   */
  public ChipSeqHMMAnalyzer(String[] args) throws IOException {

    this.args = args;
    ArgParser ap = new ArgParser(args);

    try {
      Pair<Organism, Genome> pair = Args.parseGenome(args);
      if (pair == null) {
        // Make fake genome... chr lengths provided???
        if (ap.hasKey("geninfo")) {
          genome = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
        }
        else {
          logger.fatal("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");
          printError();
          System.exit(1);
        }
      }
      else {
        genome = pair.cdr();
      }
    }
    catch (NotFoundException e) {
      e.printStackTrace();
    }
	readLength = Args.parseInteger(args,"readlen",readLength);
    //load experiments
    loadExperiments();
    
    DeepSeqExpt signal = experiments.get(0).car();
    DeepSeqExpt control = experiments.get(0).cdr();

    Double wceScaleFactor = null;
    if (control != null) {
      //TODO is there a better way to do this?
      wceScaleFactor = new Double(Math.min(signal.getWeightTotal() / control.getWeightTotal(), 1.0));      
    }


    //set up binding model
    String fragDistFilename = Args.parseString(args, FRAG_DIST_FILENAME_KEY, FRAG_DIST_FILENAME_DFLT);
    int model_min = Args.parseInteger(args, MODEL_MIN_KEY, MODEL_MIN_DFLT);
    int model_max = Args.parseInteger(args, MODEL_MAX_KEY, MODEL_MAX_DFLT);
    File fragDistFile = new File(fragDistFilename);
    if (!fragDistFile.isFile()) {
      logger.fatal("Cannot find binding model file");
      System.exit(1);
    }
    BindingModel model = new BindingModel(fragDistFile, model_min, model_max);        

    
    //read selected regions
    String regionFilename = Args.parseString(args, SELECTED_REGION_FILENAME_KEY, SELECTED_REGION_FILENAME_DFLT);
    String regionFormat = Args.parseString(args, SELECTED_REGION_FILE_FORMAT_KEY, SELECTED_REGION_FILE_FORMAT_DFLT);
    double minReads = Args.parseDouble(args, MINIMUM_READS_FOR_SELECTED_REGION_KEY, MINIMUM_READS_FOR_SELECTED_REGION_DFLT);
    
    ArrayList<Region> selectedRegions = null;    
    if (regionFormat.equals("NONE")) {
      if (signal.isFromDB()) {
        logger.fatal("Cannot estimate restrict regions from DB.");
        System.exit(1);
      }
//      selectedRegions = ChipSeqAnalysisUtils.selectEnrichedRegions(genome, model.getWidth(), model.getRange(), minReads, signal, control);
      ChipSeqAnalysisUtils.saveRestrictRegions(selectedRegions, true);
    }
    else if (regionFormat.equals("MACS")) {
      File regionFile = new File(regionFilename);
      List<MACSPeakRegion> peaks = MACSParser.parseMACSOutput(regionFile.getAbsolutePath(), genome);
      selectedRegions = ChipSeqAnalysisUtils.selectEnrichedMACSRegions(model.getRange(), peaks);
    }
    else if (regionFormat.equals("StatPeak")) {
      selectedRegions = ChipSeqAnalysisUtils.loadRegionsFromFile(genome, regionFilename, false, model.getRange());
    }
    else if (regionFormat.equals("Regions")) {
      selectedRegions = ChipSeqAnalysisUtils.loadRegionsFromFile(genome, regionFilename, true,model.getRange());
    }
    else {
      logger.fatal("Unable to load selected regions");      
      System.exit(1);
    }


    //read/initialize weight matrices
    List<WeightMatrix> motifPWMs = new ArrayList<WeightMatrix>();
    WeightMatrix bgPWM = null;
    String motifPWMName = Args.parseString(args, MOTIF_PWM_NAME, null);
    String motifPWMFilename = Args.parseString(args, MOTIF_PWM_FILENAME, null);
    if (motifPWMFilename != null) {
      WeightMatrix motifPWM = null;      
      List<WeightMatrix> pwms = PWMParser.readTRANSFACFreqMatrices(motifPWMFilename, "TRANSFAC");
      motifPWM = pwms.get(0);
      motifPWMs.add(motifPWM);
    }
    else if (motifPWMName != null) {
      //TODO figure out how to convert from log odds to frequency format
//      String motifPWMVersion = Args.parseString(args, MOTIF_PWM_VERSION, "");
//      WeightMatrix motifPWM = null;
//      try {
//        int pwmID = WeightMatrix.getWeightMatrixID(genome.getSpeciesDBID(), motifPWMName, motifPWMVersion);
//        motifPWM = WeightMatrix.getWeightMatrix(pwmID);
//        motifPWM.toFrequency();
//      }
//      catch (NotFoundException nfex) {
//        logger.fatal("Invalid Motif PWM Specification - name: " +motifPWMName + ", version: " + motifPWMVersion);
//        System.exit(1);
//      }
//      motifPWMs.add(motifPWM);
    }
    else {
      /**
       * TODO enumerate all n-mers and run EM for a single iteration for each
       * to determine the best candidates
       */
      logger.fatal("No Motif PWM specified - not yet implemented");
      System.exit(1);
    }
    
    String bgPWMName = Args.parseString(args, BG_PWM_NAME, null);
    String bgPWMFilename = Args.parseString(args, BG_PWM_FILENAME, null);
    if (bgPWMFilename != null) {      
        List<WeightMatrix> pwms = PWMParser.readTRANSFACFreqMatrices(bgPWMFilename,"TRANSFAC");
      bgPWM = pwms.get(0);
    }    
    else if (bgPWMName != null) {
      String bgPWMVersion = Args.parseString(args, BG_PWM_VERSION, "");      
      try {
        int pwmID = WeightMatrix.getWeightMatrixID(genome.getSpeciesDBID(), bgPWMName, bgPWMVersion);
        bgPWM = WeightMatrix.getWeightMatrix(pwmID);
      }
      catch (NotFoundException nfex) {
        logger.fatal("Invalid Motif PWM Specification - name: " + bgPWMName + ", version: " + bgPWMVersion);
        System.exit(1);
      }
    }
    else {
      bgPWM = this.makeDefaultBGPWM();
    }
    
    
    //set up hmm model for binding and motifs
    chipSeqModel = new BindingMotifHMM(signal, control, wceScaleFactor, selectedRegions, model, args);

    //initialize the HMM
    chipSeqModel.initializeHMM(motifPWMs, bgPWM);
  }


  /**
   * 
   */
  public void loadExperiments() {
    // Experiments : Load each condition expt:ctrl Pair
    long loadData_tic = System.currentTimeMillis();
    ArrayList<String> conditionNames = new ArrayList<String>();
    int exptHitCount = 0;
    int ctrlHitCount = 0;
    Vector<String> exptTags = new Vector<String>();
    for (String s : args) {
      if (s.contains("expt")) {
        if (!exptTags.contains(s)) {
          exptTags.add(s);
        }
      }
    }

    // each tag represents a condition
    for (String tag : exptTags) {
      String name = "";
      if (tag.startsWith("--db")) {
        name = tag.replaceFirst("--dbexpt", "");
        conditionNames.add(name);
      }
      else if (tag.startsWith("--rdb")) {        
        name = tag.replaceFirst("--rdbexpt", "");
        conditionNames.add(name);
      }
      else {
        name = tag.replaceFirst("--expt", "");
        conditionNames.add(name);
      }

      if (name.length() > 0) {
        logger.debug("Loading condition: " + name);
      }

      List<ChipSeqLocator> dbexpts = Args.parseChipSeq(args, "dbexpt" + name);
      List<ChipSeqLocator> dbctrls = Args.parseChipSeq(args, "dbctrl" + name);
      List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args, "rdbexpt" + name);
      List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args, "rdbctrl" + name);
      List<File> expts = Args.parseFileHandles(args, "expt" + name);
      List<File> ctrls = Args.parseFileHandles(args, "ctrl" + name);
      ArgParser ap = new ArgParser(args);
      boolean nonUnique = ap.hasKey("nonunique") ? true : false;
      String fileFormat = Args.parseString(args, "format", "ELAND");

      if (expts.size() > 0 && dbexpts.size() == 0 && rdbexpts.size() == 0) {
        DeepSeqExpt e = new DeepSeqExpt(genome, expts, nonUnique, fileFormat, exptHitCount);
        DeepSeqExpt c = new DeepSeqExpt(genome, ctrls, nonUnique, fileFormat, ctrlHitCount);
        experiments.add(new Pair<DeepSeqExpt, DeepSeqExpt>(e, c));
        exptHitCount += e.getHitCount();
        ctrlHitCount += c.getHitCount();
      }
      else if (dbexpts.size() > 0 && expts.size() == 0) {
        experiments.add(new Pair<DeepSeqExpt, DeepSeqExpt>(new DeepSeqExpt(genome, dbexpts, "db", readLength), 
            new DeepSeqExpt(genome, dbctrls, "db", readLength)));
      }
      else if (rdbexpts.size() > 0 && expts.size() == 0) {
        experiments.add(new Pair<DeepSeqExpt, DeepSeqExpt>(new DeepSeqExpt(genome, rdbexpts, "readdb", readLength),
            new DeepSeqExpt(genome, rdbctrls, "readdb", readLength)));
      }
      else {
        logger.fatal("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
        printError();
        System.exit(1);
      }
    }
    logger.debug("Done creating DeepSeqExpt objects");
  }


  public void runMixtureModel() {

    String peakFileName = chipSeqModel.getOutName();

    // mixModel.countNonSpecificReads();

    chipSeqModel.execute();
    chipSeqModel.setOutName(peakFileName);
    
    chipSeqModel.printFeatures(); 

    // mixModel.addAnnotations();
    //TODO Need to add this to the model or it won't do anything 
    this.printInsignificantFeatures(peakFileName); 

    // mixModel.printPeakSequences();    
    logger.debug("Finished! Binding events are printed to: " + chipSeqModel.getOutName());
  }


  public WeightMatrix makeDefaultBGPWM() {
    WeightMatrix bgPWM = new WeightMatrix(1);
    for (int i = 0; i < WeightMatrix.allLetters.length; i++) {
      bgPWM.matrix[0][WeightMatrix.allLetters[i]] = (float)0.25;
    }
    return bgPWM;
  }

  /**
   * Configure log4j
   */
	public static void configureLogging() {
		ClassLoader loader = ChipSeqHMMAnalyzer.class.getClassLoader();
		PropertyConfigurator.configure(loader.getResource("edu/mit/csail/cgs/utils/config/log4j.properties"));
	}

	
	public static String[] setTestArgs() {
	  String[] args = new String[] {
	      "--species", "Mus musculus;mm8",
	      "--geninfo", "Mus musculus;mm8",
	      "--rdbexptX", "Sing_CTCF_ES;bowtie_unique",
	      "--rdbctrlX", "Sing_GFP_ES;bowtie_unique",
	      "--region_file", "hand_test_region2.txt",
	      "--region_file_format", "Regions",
	      "--minReads", "0",
	      "--frag_file", "oct4.shear.txt",
	      "--model_min", "-200",
	      "--model_max", "500",
	      "--update_motif", "FALSE",
	      "--motif_PWM_filename", "hand_test_motif.txt",
	      "--motif_binding_freq", "1.0"
	  };
	  return args;
	}
  
  public static void main(String[] args) throws Exception {
  	ChipSeqHMMAnalyzer.configureLogging();
  	//Logger.getRootLogger().setLevel(Level.DEBUG);
    logger.debug("Starting ChipSeqHMMAnalyzer");
    //args = ChipSeqHMMAnalyzer.setTestArgs();
    ChipSeqHMMAnalyzer analyzer = new ChipSeqHMMAnalyzer(args);
    analyzer.runMixtureModel();
  }


  /**
   * Command-line help
   */
  public void printError() {
    System.err.println("Usage:\n " + "ChipSeqAnalyzer \n" + "Using with Gifford Lab DB:\n"
        + "  --species <organism name;genome version>" + "  --dbexptX <IP expt (X is condition name)> "
        + "  --dbctrlX <background expt (X is condition name)> \n" + "Using wih flat-files:\n"
        + "  --geninfo <file with chr name/length pairs> "
        + "  --exptX <aligned reads file for expt (X is condition name)> "
        + "  --ctrlX <aligned reads file for ctrl (X is condition name)> "
        + "  --format <format of above files (default ELAND)> \n"
        + "  --readlen <readlength> \n"+ "");
  }

  /*
   * command line example --property "ChipSeqAnalysis.Ctcf.properties.txt"
   * --species "Mus musculus;mm8" --dbexptCtcf
   * "Chen08_Solexa_CTCF_ES;ELAND_unique_26" --dbctrlCtcf
   * "Chen08_Solexa_GFP_ES;ELAND_unique_26" --out "CTCF_Sing_DB"
   * 
   * --property "ChipSeqAnalysis.Ctcf.properties.txt" --species
   * "Mus musculus;mm8" --exptCtcf "Sing_ES_CTCF_all.bowtie.align" --ctrlCtcf
   * "Sing_ES_GFP_all.bowtie.align" --format "BOWTIE" --out "CTCF_Sing_BOWTIE"
   */
  
  
  //snipped from BindingMixture
  private void printInsignificantFeatures(String peakFileName){
//    try{
//      String fname = peakFileName+"_Insignificant.peaks.txt"; 
//      FileWriter fw = new FileWriter(fname);
//      boolean first=true;
//      for(Feature f : chipSeqModel.getInsignificantFeatures()){
//        if(first){
//          fw.write(f.headString()); 
//          first=false;
//        }
//        fw.write(f.toString());
//      }
//      fw.close();
//    } catch (IOException e) {
//      e.printStackTrace();
//    }
  }
}
