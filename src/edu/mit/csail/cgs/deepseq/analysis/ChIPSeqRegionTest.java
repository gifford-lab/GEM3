package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class ChIPSeqRegionTest{
	private final static int MAXREAD = 1000000;
	private ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> experiments = new ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>>();
	private Genome gen;
	private ArrayList<String> conditionNames = new ArrayList<String>();
	private ArrayList<Point> coords;
	private int winSize = 1000;
	private double ratio = 1;
	
	ChIPSeqRegionTest(String[] args){
		ArgParser ap = new ArgParser(args);
		
		try {
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			if(pair==null){
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo")){
					gen = new Genome("Genome", new File(ap.getKeyValue("geninfo")), true);
	        	}else{
	        		System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs."); printError();System.exit(1);
	        	}
			}else{
				gen = pair.cdr();
			}
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		
		winSize = Args.parseInteger(args, "win", winSize);
		ratio = Args.parseDouble(args, "ratio", ratio);
		
		// load the list of coordinates
		String coord_file = Args.parseString(args, "coord", null);
		coords = CommonUtils.loadCgsPointFile(coord_file, gen);		
		
        //Experiments : Load each condition expt:ctrl Pair
        Vector<String> exptTags=new Vector<String>();
        for(String s : args)
        	if(s.contains("expt"))
        		if(!exptTags.contains(s))
        			exptTags.add(s);
    	
        // each tag represents a condition
        for(String tag : exptTags){
        	String name="";
        	if(tag.startsWith("--rdb")){
        		name = tag.replaceFirst("--rdbexpt", ""); 
        		conditionNames.add(name);
        	}

        	if(name.length()>0)
        		System.out.println("Init loading condition: "+name);
        	List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt"+name);
        	List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args,"rdbctrl"+name);
        	
        	if(rdbexpts.size()>0){
	        	experiments.add(new Pair<DeepSeqExpt,DeepSeqExpt>(new DeepSeqExpt(gen, rdbexpts, "readdb", -1),new DeepSeqExpt(gen, rdbctrls, "readdb", -1)));
	        }else{
	        	System.err.println("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");
	        	printError();
	        	System.exit(1);
	        }
        }
        

	}
    public void runTest(){
		for (int i=0;i<experiments.size();i++){
			Pair<DeepSeqExpt, DeepSeqExpt> pair = experiments.get(i);
			DeepSeqExpt ip = pair.car();
			DeepSeqExpt ctrl = pair.cdr();
			System.out.println("\nCounting Experiment "+conditionNames.get(i)+" ...");
			for (Point p: coords){
				Region r = p.expand(winSize);
				int ipCount = ip.countHits(p.expand(winSize));
				int ctrlCount = ctrl.countHits(p.expand(winSize));
				double negLog = -Math.log10(StatUtil.binomialPValue(ctrlCount, ipCount+ctrlCount, ratio));
				System.out.println(String.format("%s\t%d\t%d\t%.2f", p.toString(), ipCount, ctrlCount, negLog));
			}
			System.out.println("\nDone! \n");
		}
    }	


	public static void main(String[] args){
		ChIPSeqRegionTest test = new ChIPSeqRegionTest(args);
		test.runTest();
	}

	/**
	 * Command-line help
	 */
	public void printError() {
		System.err.println("Usage:\n " +
                "ChIPSeqRegionTest \n" +
                "Using with Gifford Lab ReadDB:\n" +
                "  --species <organism name;genome version>\n"+
                "  --rdbexptX <IP expt (X is condition name)>\n" +
                "  --rdbctrlX <background expt (X is condition name)> \n" +
            	"");		
	}


	/* command line example 
--species "Mus musculus;mm8" 
--rdbexptCtcf "Sing_CTCF_ES;bowtie_unique" 
--rdbctrlCtcf  "Sing_GFP_ES;bowtie_unique"
	 */
}
