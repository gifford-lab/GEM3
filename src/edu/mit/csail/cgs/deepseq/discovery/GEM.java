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
	public final static String GEM_VERSION = "2.5";
	private String[] args;
	private Genome genome;
    private KPPMixture mixture;
	
	public GEM(String[] args) throws NotFoundException {
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
                genome = new Genome("Genome", new File(genomeString), true);
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

		// esitmate chrom sizes for the genome if it is not provided with --g option
        if(genome==null){
    		System.out.println("Estimating chromosome sizes from all read files (skip this step by adding a --g option)...\n");
			HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
        	// each tag represents a condition
            for(String tag : exptTags){
        		if(!tag.startsWith("--rf") && (!tag.startsWith("--rdb"))){
    	        	String name = tag.replaceFirst("--expt", ""); 
    	        	List<File> expts = Args.parseFileHandles(args, "expt"+name);
    	        	List<File> ctrls = Args.parseFileHandles(args, "ctrl"+name);  
    	            String fileFormat = Args.parseString(args, "f", "BED").toUpperCase();
    	            if (fileFormat.equals("BAM"))
    	            	fileFormat = "SAM";
    	            if(expts.size()>0){
    	                DeepSeqExpt e = new DeepSeqExpt(expts, fileFormat);
    	                DeepSeqExpt c = new DeepSeqExpt(ctrls, fileFormat);
    	                
    	    			Map<String, Integer> currMap = e.getGenome().getChromLengthMap();
    	    			for(String s: currMap.keySet()){
    	    				if(!chrLenMap.containsKey(s) || chrLenMap.get(s)<currMap.get(s)+1000)
    	    					chrLenMap.put(s, currMap.get(s)+1000);
    	    			}
    	    			currMap = c.getGenome().getChromLengthMap();
    	    			for(String s: currMap.keySet()){
    	    				if(!chrLenMap.containsKey(s) || chrLenMap.get(s)<currMap.get(s))
    	    					chrLenMap.put(s, currMap.get(s));
    	    			}
    	            }
        		}
        	}
			genome=new Genome("Genome", chrLenMap);
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
	            if (fileFormat.equals("BAM"))
	            	fileFormat = "SAM";
	
	            if(expts.size()>0 && rdbexpts.size()==0){
	                int readLength = -1;	// For file, read length will be obtained from the data
	                DeepSeqExpt e = new DeepSeqExpt(genome, expts, nonUnique, fileFormat, readLength);
	                DeepSeqExpt c = new DeepSeqExpt(genome, ctrls, nonUnique, fileFormat, readLength);
	                experiments.add(new Pair<DeepSeqExpt,DeepSeqExpt>(e,c));
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
	
    public void runMixtureModel(boolean run_gem) {	    	
        
        Set<String> flags = Args.parseFlags(args);
        int GPS_round = Args.parseInteger(args,"r_gps", 2);
        int GEM_round = Args.parseInteger(args,"r_gem", 2);
        int minLeft = Args.parseInteger(args,"d_l", 300);
        int minRight = Args.parseInteger(args,"d_r", 200);
        boolean not_update_model = flags.contains("not_update_model");
        /**
         ** Simple GPS event finding without sequence information
         **/
	    int round = 0;
        String filePrefix = mixture.getOutName();
        String prefix = new File(filePrefix).getName();
        mixture.setOutName(filePrefix+"_"+round);
        System.out.println("\n============================ Round "+round+" ============================");
        
        mixture.execute();
        if (!not_update_model){
	        boolean constant_model_range = Args.parseFlags(args).contains("constant_model_range");
	        if (!constant_model_range){
	            Pair<Integer, Integer> newEnds = mixture.getModel().getNewEnds(minLeft, minRight);
	            mixture.updateBindingModel(newEnds.car(), newEnds.cdr(), filePrefix+"_"+(round+1));
	        }
	        else
	            mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax(), filePrefix+"_"+(round+1));
        }
        if (Args.parseFlags(args).contains("process_all_regions")){
	        mixture.printFeatures(round);
	        mixture.printFilteredFeatures(round);
	        mixture.printInsignificantFeatures(round);	        
        	if (Args.parseFlags(args).contains("refine_regions"))
            	mixture.refineRegions();
        }
        mixture.releaseMemory();
        
        while (round+1<GPS_round){
            round++;
            System.out.println("\n============================ Round "+round+" ============================");
            mixture.setOutName(filePrefix+"_"+round);
            mixture.execute();
            if (!not_update_model)
            	mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax(), filePrefix+"_"+(round+1));
            mixture.printFeatures(round);
            mixture.printFilteredFeatures(round);
            mixture.printInsignificantFeatures(round);
	        mixture.releaseMemory();
        }

    	File currentFolder = new File(filePrefix).getParentFile().getParentFile();
    	String path = new File(currentFolder, prefix).getAbsolutePath();
        CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.txt", path+"_GPS_events.txt");
        if (Args.parseFlags(args).contains("outNP"))
        	CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.narrowPeak", path+"_GPS_events.narrowPeak");
        
        /**
         ** GEM event finding with kmer motif positional prior (KPP)
         **/    	
        if (run_gem){
        	// initialize first set of kmers from GPS result
	        int returnValue = mixture.initKMAC();	
	        if (returnValue < 0){					// this could happen if no k value can be found to give good motif
	        	mixture.plotAllReadDistributions(mixture.allModels, mixture.outName);
	            mixture.closeLogFile();
	            
	            if (returnValue == -1)
	            	System.out.println("\nMotif can not be found!\n\nGPS analysis results are printed to:");
	            else if (returnValue == -2)
	            	System.out.println("\nBinding event can not be found!\n\nGPS analysis results are printed to:");
	            
	            System.out.println(path+"_GPS_events.txt\n"+
		        		path+"_result.htm\n" +
		        		path+"_outputs (folder with all other files)\n");
		        
		        String htmName = prefix+"_outputs/"+prefix+"_"+round+"_result.htm";
		        String html = "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN'><html><head><title>Redirect</title><meta http-equiv='REFRESH' content='0;url="+
		        	htmName+"'></HEAD><BODY>If your browser did not redirect, <a href='"+
		        	htmName+"'>click here for GPS Result</a>.</BODY></HTML>";
		        CommonUtils.writeFile(path+"_result.htm", html);
	            return;
	        }
	        	
            for (int i=1;i<GEM_round;i++){
				round++;			
	            System.out.println("\n============================ Round "+round+" ============================");
	            mixture.setOutName(filePrefix+"_"+round);
	            mixture.execute();   
	            if (!not_update_model)
	            	mixture.updateBindingModel(-mixture.getModel().getMin(), mixture.getModel().getMax(), filePrefix+"_"+(round+1));
		        mixture.printFeatures(round);
		        mixture.printFilteredFeatures(round);
		        mixture.printInsignificantFeatures(round);
		        mixture.releaseMemory();
		        mixture.runKMAC(Args.parseInteger(args,"k_win", 61));// Note: KPPMixture also has args parsing, keep default value the same
            }
            int winSize = Args.parseInteger(args,"k_win2", -1);
            if (winSize!=-1){
	            System.out.println("\n============== Finding motif for "+prefix+"_"+(round+1)+", window size="+winSize+" =============\n");
	            mixture.setOutName(filePrefix+"_"+(round+1));
		        mixture.runKMAC(winSize);	// Note: KPPMixture also has args parsing, keep default value the same
            }
        }
        
        mixture.plotAllReadDistributions(mixture.allModels, mixture.outName);
        mixture.closeLogFile();        

        if (run_gem){
	        System.out.println("\nFinished! GEM analysis results are printed to:\n"+
	        		path+"_GEM_events.txt\n"+
	        		path+"_KSM.txt\n"+
	        		path+"_PFM.txt\n"+
	        		path+"_result.htm\n" +
	        		path+"_outputs (folder with all other files)\n");
	        CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.txt", path+"_GEM_events.txt");
	        if (Args.parseFlags(args).contains("outNP"))
	        	CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.narrowPeak", path+"_GEM_events.narrowPeak");
	        CommonUtils.copyFile(filePrefix+"_"+round+"_PFM.txt", path+"_PFM.txt");
	        CommonUtils.copyFile(filePrefix+"_"+round+"_KSM.txt", path+"_KSM.txt");
	        String htmName = prefix+"_outputs/"+prefix+"_"+round+"_result.htm";
	        String html = "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN'><html><head><title>Redirect</title><meta http-equiv='REFRESH' content='0;url="+
	        	htmName+"'></HEAD><BODY>If your browser did not redirect, <a href='"+
	        	htmName+"'>click here for GEM Result</a>.</BODY></HTML>";
	        CommonUtils.writeFile(path+"_result.htm", html);
        }
        else{
            System.out.println("\nFinished! GPS analysis results are printed to:\n"+
            		path+"_GPS_events.txt\n"+
	        		path+"_outputs (folder with all other files)\n");
	        CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.txt", path+"_GPS_events.txt");
	        if (Args.parseFlags(args).contains("outNP"))
	        	CommonUtils.copyFile(filePrefix+"_"+round+"_GEM_events.narrowPeak", path+"_GPS_events.narrowPeak");
        }
    }
    	
    public static void main(String[] args) throws Exception {
        long tic = System.currentTimeMillis();
        
        System.out.println("\nWelcome to GEM (version "+GEM_VERSION+")!");
        System.out.println("\nPlease cite: \nYuchun Guo, Shaun Mahony, David K. Gifford (2012) PLoS Computational Biology 8(8): e1002638. \nHigh Resolution Genome Wide Binding Event Finding and Motif Discovery Reveals Transcription Factor Spatial Binding Constraints. \ndoi:10.1371/journal.pcbi.1002638\n");
        System.out.println("Gifford Laboratory at MIT (http://cgs.csail.mit.edu/gem/).\n");
    	
        boolean run_gem = false;
    	if (Args.parseInteger(args,"k", -1)!=-1 || Args.parseInteger(args,"k_min", -1)!=-1 || Args.parseInteger(args,"kmin", -1)!=-1
    			|| Args.parseString(args, "seed", null)!=null)
    		run_gem = true;
    	else
    		System.err.println("Warning: GEM did not see options (--k, --k_min & --k_max, or --seed) to run motif discovery. It will run GPS and stop!");

        GEM gem = new GEM(args);
        
        gem.runMixtureModel(run_gem);
        gem.close();
        
        System.out.println("\nTotal running time: "+CommonUtils.timeElapsed(tic)+"\n");
    }

    /**
     * Command-line help
     */
    public static void printHelp() {
        System.err.print("" +
                         "GEM command line options (see more options at our website)\n" +
                         "   Required parameters:\n" +
                         "      --d <read spatial distribution file>\n" +
                         "      --exptX <aligned read file for expt (X is condition name)>\n" +
                         "   Required GEM motif discovery parameters, optional for GPS-only analysis:\n" +
                         "      --k <length of the k-mer for motif finding, use --k or (--kmin & --kmax)>\n" +
                         "      --k_min <min value of k, e.g. 6>\n" +
                         "      --k_max <max value of k, e.g. 13>\n" +
                         "      --seed <exact k-mer string to jump start k-mer set motif discovery>\n" +
                         "      --genome <the path to the genome sequence directory, for motif finding>\n" +
                         "   Optional parameters:\n" +
                         "      --ctrlX <aligned reads file for ctrl (for each condition, ctrlX should match exptX)>\n" +
                         "      --g <genome chrom.sizes file with chr name/length pairs>\n" +
                         "      --f <read file format, BED/SAM/BOWTIE/ELAND/NOVO (default BED)>\n" +
                         "      --s <size of mappable genome in bp (default is estimated from genome chrom sizes)>\n" +
                         "      --a <minimum alpha value for sparse prior (default is esitmated from the whole dataset coverage)>\n" +
                         "      --q <significance level for q-value, specify as -log10(q-value) (default=2, q-value=0.01)>\n" +
                         "      --t <maximum number of threads to run GEM in paralell (default=#CPU)>\n" +
                         "      --out <output folder name and file name prefix>\n" +
                         "      --k_seqs <number of binding events to use for motif discovery (default=5000)>\n" +
                         "   Optional flags: \n" +
                         "      --fa use a fixed user-specified alpha value for all the regions\n" +
                         "      --help print this help information and exit\n" +
//                         "\n   Output format:\n" +
//                         "      The output file contains eight fields in a tab-delimited file:\n" +
//                         "        - Binding event coordinate\n" +
//                         "        - IP read count\n" +
//                         "        - Control read count\n" +
//                         "        - Fold enrichment (IP/Control)\n" +                
//                         "        - P-value\n" +
//                         "        - Q-value (multiple hypothesis corrected)\n"+
//                         "        - Shape deviation from the empirical read distribution (log10(KL))\n" +
//                         "        - Shape deviation between IP vs Control (log10(KL))\n" +
                         "\n"+"Example: java -Xmx10G -jar gem.jar --d Read_Distribution_default.txt --g mm8.chrom.sizes --genome your_path/mm8 --s 2000000000 --expt SRX000540_mES_CTCF.bed --ctrl SRX000543_mES_GFP.bed --f BED --out mouseCTCF --k_min 6 --k_max 13\n");	
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

