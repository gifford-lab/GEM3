package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class BEDFileWriter{
	private final static int MAXREAD = 1000000;
	private ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> experiments = new ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>>();
	private Genome gen;
	private double fraction = 1;
	private ArrayList<String> conditionNames = new ArrayList<String>();
	private String chr;
	
	BEDFileWriter(String[] args){
		ArgParser ap = new ArgParser(args);
		gen = CommonUtils.parseGenome(args);
		fraction = Args.parseDouble(args,"fraction",fraction);
		chr = Args.parseString(args, "chr", null);
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
    public void writeBED(){
		for (int i=0;i<experiments.size();i++){
			Pair<DeepSeqExpt, DeepSeqExpt> pair = experiments.get(i);
			DeepSeqExpt ip = pair.car();
			DeepSeqExpt ctrl = pair.cdr();
			String name_ip = conditionNames.get(i)+(chr!=null?"_chr"+chr:"")+"_ip"+ (fraction==1?"":("_"+fraction))+".bed";
			String name_ctrl = conditionNames.get(i)+(chr!=null?"_chr"+chr:"")+"_ctrl"+ (fraction==1?"":("_"+fraction))+".bed";
			// clean the target file if exists
			boolean fileCleaned = resetFile(name_ip);
			if (!fileCleaned){
				System.err.println(name_ip+" can not be reset. Skipped!");
				continue;
			}
			fileCleaned = resetFile(name_ctrl);
			if (!fileCleaned){
				System.err.println(name_ctrl+" can not be reset. Skipped!");
				continue;
			}
			
			System.out.println("\nWriting Experiment "+conditionNames.get(i)+" ...");
			List<String> chroms = new ArrayList<String>();
			if (chr!=null)
				chroms.add(chr);
			else
				chroms =gen.getChromList();
			Collections.sort(chroms);
			for (String chrom: chroms){
				System.out.println("Writing Chomosome "+chrom+" ...");
				// load  data for this chromosome.
				int length = gen.getChromLength(chrom);
				Region wholeChrom = new Region(gen, chrom, 0, length);
				int count = Math.max(ip.countHits(wholeChrom), ctrl.countHits(wholeChrom));
				ArrayList<Region> chunks = new ArrayList<Region>();
				// if there are too many reads in a chrom, read smaller chunks
				if (count>MAXREAD){
					int chunkNum = count/MAXREAD+1;
					int chunkLength = length/chunkNum;
					int start = 0;
					while (start<=length){
						int end = Math.min(length, start+chunkLength-1);
						Region r = new Region(gen, chrom, start, end);
						start = end+1;
						chunks.add(r);
					}
				}else
					chunks.add(wholeChrom);
					
				for (Region chunk: chunks){
					writeFile(name_ip,ip.getBED_StrandedReads(chunk, '+', fraction));
					writeFile(name_ip,ip.getBED_StrandedReads(chunk, '-', fraction));
					if (ctrl!=null){
						writeFile(name_ctrl, ctrl.getBED_StrandedReads(chunk, '+', fraction));
						writeFile(name_ctrl, ctrl.getBED_StrandedReads(chunk, '-', fraction));
					}
				}
			}
			System.out.println("\nDone! \n"+name_ip+"\n"+name_ctrl);
		}
    }	
	public static void writeFile(String fileName, String text){
		try{
			FileWriter fw = new FileWriter(fileName, true); //append file
			fw.write(text);
			fw.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}	
	// reset the target file if exists
	private boolean resetFile(String name){
		boolean exists = (new File(name)).exists();
		if (exists) {
			boolean success = (new File(name)).delete();
			if (!success) {
			    return false;
			}
		} 
		return true;
	}
	
	public static void main(String[] args){
		BEDFileWriter writer = new BEDFileWriter(args);
		writer.writeBED();
	}

	/**
	 * Command-line help
	 */
	public void printError() {
		System.err.println("Usage:\n " +
                "BEDFileWriter \n" +
                "Using with Gifford Lab Read DB:\n" +
                "  --species <organism name;genome version>\n"+
                "  --rdbexptX <IP expt (X is condition name)>\n" +
                "  --rdbctrlX <background expt (X is condition name)> \n" +
                "Optional\n" +
                "  --fraction <fraction of reads to output>\n" +
            	"");		
	}


	/* command line example 
--species "Mus musculus;mm8" 
--rdbexptCtcf "Sing_CTCF_ES;bowtie_unique" 
--rdbctrlCtcf  "Sing_GFP_ES;bowtie_unique"
	 */
}
