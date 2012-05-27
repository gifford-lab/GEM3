package edu.mit.csail.cgs.deepseq.discovery;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;


public class GEM {
	public final static String GEM_VERSION = "0.9";
	private String[] args;
	private Genome genome;
    private KPPMixture mixture;
	
	public GEM(String[] args) throws NotFoundException {
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
            String genomeString = Args.parseString(args,"g",null);		// text file with chrom lengths
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
        	if(tag.startsWith("--rf")){
        		name = tag.replaceFirst("--rfexpt", ""); 
        		conditionNames.add(name);
        	}
        	else {
        		System.out.println("Loading data...");
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
	        	boolean nonUnique = flags.contains("nonunique");
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
	                System.out.println("    done: "+CommonUtils.timeElapsed(loadData_tic));
	            }
	        }
        }
        
        try{
            mixture = new KPPMixture(genome, experiments, conditionNames, args);
        } catch(Exception ex){
            for(Pair<DeepSeqExpt,DeepSeqExpt> e : experiments){
                e.car().closeLoaders();
                e.cdr().closeLoaders();
            }
            ex.printStackTrace();
        }
    } //end of GPS constructor
	
    public void runMixtureModel() {		
    	boolean run_gem = false;
    	if (Args.parseInteger(args,"k", -1)!=-1 || Args.parseInteger(args,"k_min", -1)!=-1
    			|| Args.parseString(args, "seed", null)!=null)
    		run_gem = true;
        double kl=10;
		
//        run_gem = false;		// DO NOT RUN GEM, for GPS v1.1 release
        Set<String> flags = Args.parseFlags(args);
        int GPS_round = Args.parseInteger(args,"r_gps", 2);
        int GEM_round = Args.parseInteger(args,"r_gem", 2);
        int minLeft = Args.parseInteger(args,"d_l", 300);
        int minRight = Args.parseInteger(args,"d_r", 200);
        /**
         ** Simple GPS1 event finding without sequence information
         **/
	    int round = 0;
        String filePrefix = mixture.getOutName();
        String prefix = new File(filePrefix).getName();
        mixture.setOutName(filePrefix+"_"+round);
        System.out.println("\n============================ Round "+round+" ============================");
        mixture.execute();
        mixture.printFeatures(round);
        mixture.printFilteredFeatures(round);
        mixture.printInsignificantFeatures(round);
        mixture.refineRegions();
        
        while (round+1<GPS_round){
            round++;
            System.out.println("\n============================ Round "+round+" ============================");
            mixture.setOutName(filePrefix+"_"+round);
            if (round==1){            	
                boolean noChange = Args.parseFlags(args).contains("constant_model_range");
                if (!noChange){
                    Pair<Integer, Integer> newEnds = mixture.getModel().getNewEnds(minLeft, minRight);
                    kl = mixture.updateBindingModel(newEnds.car(), newEnds.cdr());
                }
                else
                    kl = mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax());
            }
            else
                kl = mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax());
            mixture.execute();
            mixture.printFeatures(round);
            mixture.printFilteredFeatures(round);
            mixture.printInsignificantFeatures(round);
        }
        
        /**
         ** GPS2 event finding with kmer positional prior (KPP)
         **/   
        if (run_gem){
        	// initialize first set of kmers from GPS result
	        int returnValue = mixture.initKMF();	
	        if (returnValue == -1){					// this could happen if no k value can be found to give good motif
	        	mixture.plotAllReadDistributions();
	            mixture.closeLogFile();
	            return;
	        }
	        	
            for (int i=1;i<GEM_round;i++){
				round++;			
	            System.out.println("\n============================ Round "+round+" ============================");
	            mixture.setOutName(filePrefix+"_"+round);
	            mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax());
	            mixture.execute();   
		        mixture.printFeatures(round);
		        mixture.printFilteredFeatures(round);
		        mixture.printInsignificantFeatures(round);
		        mixture.runKMF(Args.parseInteger(args,"k_win", 61));// Note: KPPMixture also has args parsing, keep default value the same
            }
            int winSize = Args.parseInteger(args,"k_win2", -1);
            if (winSize!=-1){
	            System.out.println("\n============== Finding motif for "+prefix+"_"+(round+1)+", large window size="+winSize+" =============\n");
	            mixture.setOutName(filePrefix+"_"+(round+1));
		        mixture.runKMF(winSize);	// Note: KPPMixture also has args parsing, keep default value the same
            }
        }
        
        mixture.plotAllReadDistributions();
        mixture.closeLogFile();
        
//        round --;

    	File currentFolder = new File(filePrefix).getParentFile().getParentFile();
    	String path = new File(currentFolder, prefix).getAbsolutePath();
        if (run_gem){
	        System.out.println("\nFinished! GEM analysis results are printed to:\n"+
	        		path+"_GEM_events.txt\n"+
	        		path+"_result.htm\n" +
	        		path+"_outputs (folder)\n");
	        CommonUtils.copyFile(filePrefix+"_"+GPS_round+"_GEM_events.txt", path+"_GEM_events.txt");
	        String htmName = prefix+"_outputs/"+prefix+"_"+GPS_round+"_result.htm";
	        String html = "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN'><html><head><title>Redirect</title><meta http-equiv='REFRESH' content='0;url="+
	        	htmName+"'></HEAD><BODY>If your browser did not redirect, <a href='"+
	        	htmName+"'>click here for GEM Result</a>.</BODY></HTML>";
	        CommonUtils.writeFile(path+"_result.htm", html);
        }
        else{
            System.out.println("\nFinished! GPS analysis results are printed to:\n"+
            		path+"_GPS_events.txt\n"+
	        		path+"_outputs (folder)\n");
	        CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.txt", path+"_GPS_events.txt");
        }
    }
    	
    public static void main(String[] args) throws Exception {
        long tic = System.currentTimeMillis();
        System.out.println("\nWelcome to GEM (version "+GEM_VERSION+")!");
        System.out.println("Gifford Laboratory at MIT (http://cgs.csail.mit.edu/gem/).\n");
        GEM gps = new GEM(args);
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
                         "      --exptX <aligned reads file for expt (X is condition name)>\n" +
                         "      --ctrlX <aligned reads file for ctrl (X is condition name)>\n" +
                         "   Optional parameters:\n" +
                         "      --f <read file format, BED/BOWTIE/ELAND/NOVO (default BED)>\n" +
                         "      --g <genome info file with chr name/length pairs>\n" +
                         "      --s <size of mappable genome in bp (default is estimated from genome info)>\n" +
                         "      --r <max rounds to refine read distribution (default=3)>\n" +
                         "      --a <minimum alpha value for sparse prior (default=6)>\n" +
                         "      --q <significance level for q-value, specify as -log10(q-value), (default=2, q-value=0.01)>\n" +
                         "      --t <maximum number of threads to run GPS in paralell, (default=#CPU)>\n" +
                         "      --out <output file base name>\n" +
                         "   Optional flags: \n" +
                         "      --fa use a fixed user-specified alpha value for all the regions\n" +
                         "      --help print help information and exit\n" +
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

