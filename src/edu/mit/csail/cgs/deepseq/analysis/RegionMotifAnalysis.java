package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class RegionMotifAnalysis {
	private Genome genome;
	private Organism org;
	private String[] args;
	private String motifString;
	private WeightMatrix motif = null;
	protected String outName="out";
	private double motifThreshold=0;	
	private int windowSize;
	private String filePath;
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		RegionMotifAnalysis analysis = new RegionMotifAnalysis(args);
		analysis.printMotifs();
	}
	
	RegionMotifAnalysis(String[] args){
	    ArgParser ap = new ArgParser(args);
	    try {
	      Pair<Organism, Genome> pair = Args.parseGenome(args);
	      if(pair==null){
	        //Make fake genome... chr lengths provided???
	        if(ap.hasKey("g")){
	          genome = new Genome("Genome", new File(ap.getKeyValue("g")), true);
	            }else{
	              System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");;System.exit(1);
	            }
	      }else{
	        genome = pair.cdr();
	        org = pair.car();
	      }
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }

		this.args = args;
		filePath = Args.parseString(args,"Regions","");
		windowSize = Args.parseInteger(args, "windowSize", 0);
		motifThreshold = Args.parseDouble(args, "motifThreshold", 10);
		outName = Args.parseString(args,"out",outName);
		// load motif
		try {
			motifString = Args.parseString(args, "motif", null);
			String motifVersion = Args.parseString(args, "version", null);
			int wmid = WeightMatrix.getWeightMatrixID(org.getDBID(), motifString, motifVersion);
			motif = WeightMatrix.getWeightMatrix(wmid);
		} 
		catch (NotFoundException e) {
			e.printStackTrace();
		}		
	}

	private void printMotifs(){
		long tic = System.currentTimeMillis();
		StringBuilder output = new StringBuilder();
		ArrayList<Region> regions = loadCgsRegionFile(genome, filePath, windowSize);
		
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		int count = regions.size();
		System.out.print("\nGetting motifs from "+count+" regions ...\n");
		for(int i=0; i<count; i++){
			if (i % 1000==0)
				System.out.println(i);
			Region r= regions.get(i);
			if (r.getWidth()< motif.length()){
				r=r.expand(motif.length(), motif.length());
			}
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			//search for whole region
			boolean motifFound=false;
			double score=0;
			for(int z=0; z<r.getWidth(); z++){		
				score = profiler.getMaxScore(z);
				if(score >= motifThreshold){
					motifFound=true;
					break;
				}
			}
			if (motifFound)
				output.append(r.getMidpoint().toString()).append(String.format("\t%.1f\n", score));
		}
		System.out.println(CommonUtils.timeElapsed(tic));
		CommonUtils.writeFile(outName, output.toString());
	}	

	// load text file in CGS Region format
	// chr:coord, e.g. 1:234234
	public static ArrayList<Region> loadCgsRegionFile(Genome g, String filename, int windowSize) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Region> rs = new ArrayList<Region>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				Region r = Region.fromString(g, line);
				if (windowSize!=0){
					r=r.getMidpoint().expand(windowSize);
				}
				rs.add(r);
			}
		}
		catch(IOException ioex) {
			//logger.error("Error parsing file", ioex);
		}
		finally {
			try {
				if (bin != null) {
					bin.close();
				}
			}
			catch(IOException ioex2) {
				//nothing left to do here, just log the error
				//logger.error("Error closing buffered reader", ioex2);
			}			
		}
		return rs;
	}

}
