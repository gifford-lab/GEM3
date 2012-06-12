package edu.mit.csail.cgs.deepseq.discovery;

import java.io.File;
import java.util.List;

import org.apache.log4j.Logger;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;

/**
 * Superclass for all single-condition feature discovery classes.
 * 
 * @author shaun
 * 
 */
public abstract class SingleConditionFeatureFinder extends FeatureFinder {

  private static final Logger logger = Logger.getLogger(SingleConditionFeatureFinder.class);
  
  protected DeepSeqExpt signal=null;

  protected DeepSeqExpt control=null;

  protected boolean noControl = true;


  // Constructors
  public SingleConditionFeatureFinder(DeepSeqExpt s) {
    this(s, null);
  }


  public SingleConditionFeatureFinder(DeepSeqExpt s, DeepSeqExpt c) {
    super(s.getGenome());
    signal = s;
    control = c;
    if (control != null) {
      noControl = false;
    }    
  }


  /**
   * This constructor gets used by StatisticalPeakFinder, but the code is probably
   * too specific to be used by any other class subclassing from this one
   */
  public SingleConditionFeatureFinder(String[] args) {
    super(args);

    ArgParser ap = new ArgParser(args);
    // Load the experiments
    List<ChipSeqLocator> dbexpts = Args.parseChipSeq(args, "dbexpt");
    List<ChipSeqLocator> dbctrls = Args.parseChipSeq(args, "dbctrl");
        List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt");
        List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args,"rdbctrl");
    List<File> expts = Args.parseFileHandles(args, "expt");
    List<File> ctrls = Args.parseFileHandles(args, "ctrl");
    boolean nonUnique = ap.hasKey("nonunique") ? true : false;
    boolean sigPairedEnd = ap.hasKey("sigpaired");
    boolean ctrlPairedEnd = ap.hasKey("ctrlpaired");
    String fileFormat = Args.parseString(args, "format", "ELAND");
        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
      signal = new DeepSeqExpt(gen, expts, nonUnique, fileFormat, (int)readLength);
    }
    else if (dbexpts.size() > 0 && expts.size() == 0) {
      signal = new DeepSeqExpt(gen, dbexpts, "db", (int)readLength);
      dbconnected = true;
    }
    else if (rdbexpts.size()>0 && expts.size() == 0){
        	signal = new DeepSeqExpt(gen, rdbexpts, "readdb", (int)readLength);
        	dbconnected=true;
    }
    else {
      logger.error("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
      printError();
      System.exit(1);
    }
   	signal.setPairedEnd(sigPairedEnd);
    if (ctrls.size() > 0 && dbctrls.size() == 0 && rdbctrls.size() == 0) {
      control = new DeepSeqExpt(gen, ctrls, nonUnique, fileFormat, (int)readLength);
      noControl = false;
      control.setPairedEnd(ctrlPairedEnd);
    }
    else if (dbctrls.size() > 0 && ctrls.size() == 0) {
      control = new DeepSeqExpt(gen, dbctrls, "db", (int)readLength);
      noControl = false;
      dbconnected = true;
      control.setPairedEnd(ctrlPairedEnd);
    } 
    else if (rdbctrls.size()>0 && ctrls.size() == 0) {
      control = new DeepSeqExpt(gen, rdbctrls, "readdb", (int)readLength); 
      noControl=false;
      dbconnected=true;
  	  control.setPairedEnd(ctrlPairedEnd);
    } 
    else {
      if (dbctrls.size() > 0 && ctrls.size() > 0) {
        logger.error("Cannot mix files and db loading yet...");
        printError();
        System.exit(1);
      }
      else {
        noControl = true;
        control = null;
      }
    }


    // Print some info
    logger.info("Signal hit count: " + (int) signal.getHitCount() + ", weight: " + (int) signal.getWeightTotal());
    if (!noControl) {
      logger.info("Control hit count: " + (int) control.getHitCount() + ", weight: "
          + (int) control.getWeightTotal());
    }
    logger.info("Genome size: " + genomeLen);

  }

  public DeepSeqExpt getSignal()  { return signal; }
  
  public DeepSeqExpt getControl() { return control; }

  //Call this method before exiting
  public void cleanup(){
	 if(signal!=null)
		 signal.closeLoaders();
	 if(control!=null)
		 control.closeLoaders();
  }
}
