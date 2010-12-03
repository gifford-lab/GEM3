package edu.mit.csail.cgs.deepseq.discovery;

import java.io.*;
import java.util.*;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;

import cern.jet.random.Poisson;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.*;
import edu.mit.csail.cgs.deepseq.multicond.MultiIndependentMixtureCounts;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadCache;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Utils;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.models.data.DataRegression;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;


public class GPS {
	public final static String GPS_VERSION = "1.0";
	private String[] args;
	private Genome genome;
    private GPSMixture mixture;
	
	public GPS(String[] args) throws NotFoundException {
        init();
        parseArgs(args);
    }
    public void init() {        
    }
    public void parseArgs(String[] args) throws NotFoundException {
		this.args = args;
		Set<String> flags = Args.parseFlags(args);
        if (flags.contains("help")) {
			printHelp();
			System.exit(1);
        }
        Pair<Organism, Genome> pair = Args.parseGenome(args);
        if(pair != null) {
            genome = pair.cdr();
        } else {
            String genomeString = Args.parseString(args,"g",null);
            if(genomeString != null){
                genome = new Genome("Genome", new File(genomeString));
            } else{
                genome=null;
            }
        }
		String modelFile = Args.parseString(args, "d", null);	// read distribution file
		if (modelFile == null){
			System.err.println("The read distribution file is required. Use --d option.\n");
			printError();
			System.exit(1);
        } else {
			File pFile = new File(modelFile);
			if(!pFile.isFile()){
				System.err.println("\nCannot find read distribution file!");
				System.exit(1);
			}
		}
		System.out.println("Loading data...");
        //Experiments : Load each condition expt:ctrl Pair
		ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> experiments = new ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>>();
		long loadData_tic = System.currentTimeMillis();
		ArrayList<String> conditionNames = new ArrayList<String>();
		int exptHitCount=0;
		int ctrlHitCount=0;
		Vector<String> exptTags=new Vector<String>();
		for(String s : args)
        	if(s.contains("expt"))
        		if(!exptTags.contains(s))
        			exptTags.add(s);
		
		if(exptTags.size()==0){
		    System.err.println("Error: No signal experiments provided.\nUse the --expt option.");
		    printError();
		    System.exit(1);
		}
        // each tag represents a condition
        for(String tag : exptTags){
        	String name="";
        	if(tag.startsWith("--rdb")){
        		name = tag.replaceFirst("--rdbexpt", ""); 
        		conditionNames.add(name);
        	}else{
        		name = tag.replaceFirst("--expt", ""); 
        		conditionNames.add(name);
        	}

        	if(name.length()>0)
        		System.out.println("    loading condition: "+name);
        	
        	List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt"+name);
        	List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args,"rdbctrl"+name);
        	List<File> expts = Args.parseFileHandles(args, "expt"+name);
        	List<File> ctrls = Args.parseFileHandles(args, "ctrl"+name);  
        	boolean nonUnique = flags.contains("nonunique") ? true : false;
            String fileFormat = Args.parseString(args, "f", "BED").toUpperCase();

            if(expts.size()>0 && rdbexpts.size()==0){
                int readLength = -1;	// For file, read length will be obtained from the data
                DeepSeqExpt e = new DeepSeqExpt(genome, expts, nonUnique, fileFormat, readLength);
                DeepSeqExpt c = new DeepSeqExpt(genome, ctrls, nonUnique, fileFormat, readLength);
                if(genome==null){
                    genome = DeepSeqExpt.combineFakeGenomes(e,c);
                    e.setGenome(genome);
                    c.setGenome(genome);
                }
                experiments.add(new Pair<DeepSeqExpt,DeepSeqExpt>(e,c));
                exptHitCount+=e.getHitCount();
                ctrlHitCount+=c.getHitCount();
            } else if(rdbexpts.size()>0 && expts.size() == 0){
                if(genome==null){
                    System.err.println("Error: the genome must be defined in order to use the Gifford Lab DB."); 
                    System.exit(1);
                }
                int readLength = -1;
                experiments.add(new Pair<DeepSeqExpt,DeepSeqExpt>(new DeepSeqExpt(genome, rdbexpts, "readdb", readLength),new DeepSeqExpt(genome, rdbctrls, "readdb", readLength)));
            } else{
                System.err.println("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
                printError();
                System.exit(1);
            }
        }
        System.out.println("    done: "+CommonUtils.timeElapsed(loadData_tic));
        try{
            mixture = new GPSMixture(genome, experiments, conditionNames, args);
        } catch(Exception ex){
            for(Pair<DeepSeqExpt,DeepSeqExpt> e : experiments){
                e.car().closeLoaders();
                e.cdr().closeLoaders();
            }
            ex.printStackTrace();
        }
    } //end of GPS constructor
	
    public void runMixtureModel() {		
        double kl=10;
        int round = 0;
        String peakFileName = mixture.getOutName();
        mixture.setOutName(peakFileName+"_"+round);
		
        int update_model_round = Args.parseInteger(args,"r", 3);
        while (kl>-5 && round<=update_model_round){
            System.out.println("\n============================ Round "+round+" ============================");
            mixture.execute();
            mixture.printFeatures();
            mixture.printFilteredFeatures();
            mixture.printInsignificantFeatures();
			
            round++;
            mixture.setOutName(peakFileName+"_"+round);

            if (round==1){
                boolean fixModelRange = Args.parseFlags(args).contains("fix_model_range");
                if (!fixModelRange){
                    Pair<Integer, Integer> newEnds = mixture.getModel().getNewEnds();
                    kl = mixture.updateBindingModel(newEnds.car(), newEnds.cdr());
                }
                else
                    kl = mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax());
            }
            else
                kl = mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax());
        }
        round--;
        mixture.setOutName(peakFileName+"_"+round);
        mixture.printFeatures();
        mixture.printInsignificantFeatures();
        mixture.printFilteredFeatures();
        mixture.plotAllReadDistributions();
        mixture.closeLogFile();
        System.out.println("Finished! Binding events are printed to: "+mixture.getOutName()+"_GPS_significant.txt");
    }
	
    public static void main(String[] args) throws Exception {
        long tic = System.currentTimeMillis();
        System.out.println("Welcome to GPS (version "+GPS_VERSION+")!");
        GPS gps = new GPS(args);
        gps.runMixtureModel();
        gps.close();
        System.out.println("\nTotal running time: "+CommonUtils.timeElapsed(tic)+"\n");
    }

    /**
     * Command-line help
     */
    public static void printHelp() {
        System.err.print("" +
                         "GPS Usage                      (more at http://cgs.csail.mit.edu/gps/)\n" +
                         //                "   Using with Gifford Lab DB:\n" +
                         //                "      --species <organism name;genome version>\n"+
                         //                "      --dbexptX <IP expt (X is condition name)>\n" +
                         //                "      --dbctrlX <background expt (X is condition name)>\n" +
                         //                "      --readlen <read length>\n" +
                         "   Required parameters:\n" +
                         "      --d <read distribution file>\n" +
                         "      --s <size of mappable genome in bp>\n" +
                         "      --exptX <aligned reads file for expt (X is condition name)>\n" +
                         "      --ctrlX <aligned reads file for ctrl (X is condition name)>\n" +
                         "   Optional parameters:\n" +
                         "      --f <read file format BED/BOWTIE/ELAND/NOVO (default BED)>\n" +
                         "      --g <genome info file with chr name/length pairs>\n" +
                         "      --r <max times to refine read distribution (default=3)>\n" +
                         "      --a <minimum alpha value for sparse prior (default=6)>\n" +
                         "      --q <significance level for q-value, specify as -log10(q-value), (default=2, q-value=0.01)>\n" +
                         "      --out <output file base name>\n" +
                         "   Optional flags: \n" +
                         "      --fa <use a fixed user-specified alpha value for all the regions>\n" +
                         "      --help <print help information and exit>\n" +
                         "\n   Output format:\n" +
                         "      The output file contains eight fields in a tab-delimited file:\n" +
                         "        - Binding event coordinate\n" +
                         "        - IP read count\n" +
                         "        - Control read count\n" +
                         "        - Fold enrichment (IP/Control)\n" +                
                         "        - P-value\n" +
                         "        - Q-value (multiple hypothesis corrected)\n"+
                         "        - Shape deviation from the empirical read distribution (log10(KL))\n" +
                         "        - Shape deviation between IP vs Control (log10(KL))\n" +
                         "\n");	
    }
    public void printError() {
        printHelp();
        StringBuffer sb = new StringBuffer();
        sb.append("\nYour input options are:\n");
        for (String arg:args){
            if (arg.trim().indexOf(" ")!=-1)
                sb.append("\"").append(arg).append("\" ");
            else
                sb.append(arg).append(" ");
        }
        System.err.println(sb.toString()+"\n");
    }

    //Cleanup the loaders
    //Call this before exiting
    public void close(){
        mixture.cleanup();
    }
	
}

