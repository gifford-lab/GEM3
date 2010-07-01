package edu.mit.csail.cgs.deepseq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.analysis.KnnAnalysis.KnnPoint;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.MACSPeakRegion;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.mit.csail.cgs.ewok.verbs.motifs.WeightMatrixScorer;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class MethodComparisonMotifAnalysis {
	private final int NOHIT_OFFSET = 999;
	
	private boolean isPreSorted = false;	
	private int windowSize = 50;	
	private int rank = 0;
	private List<Region> restrictRegions;
	
	private Genome genome;
	private Organism org;
	private String[] args;
	private String motifString;
	private WeightMatrix motif = null;
	private String outName="out";
	
	// each element in the list is for one ChIP-Seq method
	private ArrayList<String> methodNames = new ArrayList<String>();
	private ArrayList<ArrayList<Point>> peaks = new ArrayList<ArrayList<Point>>();
	private ArrayList<HashMap<Point, MotifHit>> maps = new ArrayList<HashMap<Point, MotifHit>>();	
//	ArrayList<HashMap<Point, ArrayList<MotifHit>>> maps_cisgenome = null;
	
	private ArrayList<Point> peaks_gps = new ArrayList<Point>();
	private ArrayList<Point> peaks_macs = new ArrayList<Point>();
	private ArrayList<Point> peaks_cisgenome = new ArrayList<Point>();
	HashMap<Point, ArrayList<MotifHit>> maps_gps = null;
	HashMap<Point, ArrayList<MotifHit>> maps_macs = null;
	HashMap<Point, ArrayList<MotifHit>> maps_cisgenome = null;
	HashMap<Point, MotifHit> map_gps = null;
	HashMap<Point, MotifHit> map_macs = null;
	HashMap<Point, MotifHit> map_cisgenome = null;
	
	// motif score cutoffs
	//	double[] motifThresholds = new double[]{0,4,8,12,16,20,24,28};
	double[] motifThresholds = new double[]{5.587944030761719, 7.559912204742432, 
		11.52034854888916, 12.886543273925781,	15.858697891235352}; //Ctcf mouse
	private double motifThreshold=0;	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MethodComparisonMotifAnalysis analysis = new MethodComparisonMotifAnalysis(args);
		// peak based analysis
		// We only print out the motif (at diff threshold) distances for each peak
		// spatialResolution or motifOccurence are analysis in Matlab 
//		analysis.printMotifOffsetDistance();
		
		// motif based analysis
//		analysis.motifBasedAnalysis();
		int type = Args.parseInteger(args, "analysisType", 0);
		switch(type){
		case 0:analysis.printMotifOffsets();break;
		case 1:analysis.proximalEventAnalysis();break;
		case 2:analysis.getUnaryEventList();break;
		case 3:analysis.singleMotifEventAnalysis();break;
		default: System.err.println("Unrecognize analysis type: "+type);
		}
	}
	
	MethodComparisonMotifAnalysis(String[] args){
		this.args = args;
		
	    ArgParser ap = new ArgParser(args);
	    try {
	      Pair<Organism, Genome> pair = Args.parseGenome(args);
	      if(pair==null){
	        //Make fake genome... chr lengths provided???
	        if(ap.hasKey("geninfo")){
	          genome = new Genome("Genome", new File(ap.getKeyValue("geninfo")));
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
	    	    
		// some parameters
		windowSize = Args.parseInteger(args, "windowSize", 50);
		isPreSorted = Args.parseInteger(args, "isPreSorted", 0)==1;
		rank = Args.parseInteger(args, "rank", 0);
		motifThreshold = Args.parseDouble(args, "motifThreshold", 10);
		outName = Args.parseString(args,"out",outName);
		try{
			restrictRegions = Args.parseRegions(args);
		}
		catch (NotFoundException e){
			// do nothing
		}
		
		// load motif
		try {
			motifString = Args.parseString(args, "motif", null);
			String motifVersion = Args.parseString(args, "version", null);
			int motif_species_id = Args.parseInteger(args, "motif_species_id", -1);
//			Organism org_mouse = new Organism("Mus musculus");
			int wmid = WeightMatrix.getWeightMatrixID(motif_species_id!=-1?motif_species_id:org.getDBID(), motifString, motifVersion);
			motif = WeightMatrix.getWeightMatrix(wmid);
		} 
		catch (NotFoundException e) {
			e.printStackTrace();
		}		
	}
	// This method print three text files
	// 1. sharedMotifOffset, for spatial resolution histogram
	// Only consider the motif positions that are shared by all methods, 
	// so that we are comparing at the same binding event, each row is matched
	// The output can be further processed by Matlab to get histogram of spatial resolution
	// for only one specified threshold
	
	// 2. allMotifOffset, for sensitivity and specificity analysis
	// Consider all the motif positions that are covered by at least one method, 
	// each row is matched, if not covered in a method, rank=-1, offset=999
	// The output can be further processed by Matlab 
	// for only one specified threshold
	
	// 3. rankedMotifOffset, both for motif coverage and occurrence curve
	// Each column is ranked as its original method, so each row is not matched
	private void printMotifOffsets(){
		long tic = System.currentTimeMillis();
		readPeakLists();
		
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		for (int i=0;i<methodNames.size();i++){
			minCount = Math.min(minCount, peaks.get(i).size());
			maxCount = Math.max(maxCount, peaks.get(i).size());
		}
		if (rank==0){
			rank = minCount;
		}
		// print some information
		String msg = String.format("Motif score cutoff=%.2f, window size=%d, top %d peaks.\n",
				motifThreshold, windowSize, rank);
		System.out.print(msg);
		System.out.println("\nNumber of peaks:");
		for (int i=0;i<methodNames.size();i++){
			System.out.println(methodNames.get(i)+"\t"+peaks.get(i).size());
		}
		
		// get all regions that covered by any of the peak callers, 
		// so that the motif search will be done only once
		// process by chrom to be more efficient in sorting
		ArrayList<Region> allRegions = new ArrayList<Region>();
		for (String chrom: genome.getChromList()){
			ArrayList<Region> byChrom = new ArrayList<Region>();
			for (ArrayList<Point> ps:peaks){
				for (Point p: ps){
					if (p.getChrom().equalsIgnoreCase(chrom))
						byChrom.add(p.expand(windowSize));
				}
			}
			allRegions.addAll( mergeRegions(byChrom));
		}
		System.out.println(BindingMixture.timeElapsed(tic));
				
		// get all the strong motifs in these regions
		Pair<ArrayList<Point>, ArrayList<Double>> results = getAllMotifs(allRegions, motifThreshold);
		ArrayList<Point> allMotifs = results.car();
		ArrayList<Double> allMotifScores = results.cdr();
		System.out.println(BindingMixture.timeElapsed(tic));
		System.out.printf("%n%d motifs (in %d regions).%n%n", allMotifs.size(), allRegions.size());
		
		// Get the set of motif matches for all peak calls
		
		System.out.printf("Peaks within a %d bp window to a motif:%n", windowSize);
		for (int i=0;i<peaks.size();i++){
			HashMap<Point, MotifHit> motifs = getNearestMotif2(peaks.get(i), allMotifs, allMotifScores);
			maps.add(motifs);
			System.out.printf("%s:%d\t", methodNames.get(i), motifs.keySet().size());
		}
		System.out.println();
		System.out.println(BindingMixture.timeElapsed(tic));
		
		System.out.printf("%nMotifs close to a peak:%n");
		for(int i = 0; i < methodNames.size(); i++)
			System.out.printf("%s\t%d/%d (%.1f%%)%n", methodNames.get(i), maps.get(i).size(), allMotifs.size(), 100.0*((double)maps.get(i).size())/(double)allMotifs.size());
				
		// get the intersection (motifs shared by all methods, upto top rank)
		System.out.print("\nRunning Spatial Resolution analysis ...\n");
		SetTools<Point> setTools = new SetTools<Point>();
		Set<Point> motifs_shared =  getHighRankSites(0);
		for (int i=1;i<maps.size();i++){
			Set<Point> motifs_new =  getHighRankSites(i);
			motifs_shared = setTools.intersection(motifs_shared,motifs_new);
		}
		msg = String.format("Motif score cutoff=%.2f, window size=%d, %d shared motifs in top %d peaks.",
				motifThreshold, windowSize, motifs_shared.size(), rank);
		System.out.println(msg);
		
		StringBuilder args_sb = new StringBuilder();
		for (String arg:args){
			args_sb.append(arg).append(" ");
		}
		String args_str = args_sb.toString();
		
		// output results, the spatial resolution (offset) 
		StringBuilder sb = new StringBuilder();
		sb.append(args_str+"\t"+msg+"\n");
		sb.append("MotifHit\tChrom\t");
		for (int i=0;i<methodNames.size();i++){
			sb.append(methodNames.get(i)+"\t");
		}
		sb.append("\n");
		for(Point motif:motifs_shared){
			sb.append(motif.toString()+"\t");
			sb.append(motif.getChrom()+"\t");
			for (int i=0;i<maps.size();i++){
				sb.append(maps.get(i).get(motif).offset);
				if (i==maps.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_sharedMotifOffsets_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
		
		// ************************************************************
		
		// get the union (motifs at least by one of all methods)
		Set<Point> motifs_union = maps.get(0).keySet();
		for (int i=1;i<maps.size();i++){
			Set<Point> motifs_new = maps.get(i).keySet();
			motifs_union = setTools.union(motifs_union,motifs_new);
		}
		msg = String.format("Total %d motifs from all methods.", motifs_union.size());
		System.out.println(msg);
		
		// output results, the spatial resolution (offset) 
		sb = new StringBuilder();
		sb.append(args_str+"\t"+msg+"\n");
		sb.append("MotifHit\tChrom\t");
		for (int i=0;i<methodNames.size();i++){
			sb.append(methodNames.get(i)+"_offset\t");
			sb.append(methodNames.get(i)+"_rank\t");
		}
		sb.append("\n");
		for(Point motif:motifs_union){
			sb.append(motif.toString()+"\t");
			sb.append(motif.getChrom()+"\t");
			for (int i=0;i<maps.size();i++){
				HashMap<Point, MotifHit> m = maps.get(i);
				if (m.containsKey(motif)){
					sb.append(m.get(motif).offset).append("\t");
					sb.append(m.get(motif).rank);
				}
				else{
					sb.append(NOHIT_OFFSET).append("\t");
					sb.append(-1);
				}
				if (i==maps.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_allMotifOffsets_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
		
		System.out.println();
		System.out.printf("Top peaks close to a motif:%n");
		// ************************************************************
		// motif occurrence / motif coverage (sensitivity)
		// get the motif offset of the peaks (rank indexed by each method)
		// if no motif hit, offset=NOHIT=999
		// print out the offset, for Matlab processing and plotting
		ArrayList<int[]> trueHits = new ArrayList<int[]>();
		for (int i=0; i<methodNames.size();i++){
			int[] hits = new int[maxCount];
			for (int j=0;j<maxCount;j++){
				hits[j]=NOHIT_OFFSET;
			}
			int numExistingTopHits = 0; 
			for (MotifHit hit: maps.get(i).values()){
				if (hit.rank < rank) {
					hits[hit.rank]=hit.offset;
					numExistingTopHits++;
				}
			}
			trueHits.add(hits);
			System.out.printf("%s\t%d/%d (%.1f%%)%n", methodNames.get(i), numExistingTopHits, rank, 100.0*((double)numExistingTopHits)/(double)rank);
		}		
		// output results
		sb = new StringBuilder();
		sb.append(args_str+"\t"+msg+"\n");
		sb.append("Rank\t");
		for (int i=0;i<methodNames.size();i++){
			sb.append(methodNames.get(i)+"\t");
		}
		sb.append("\n");
		for(int j=0;j<maxCount;j++){
			sb.append(j+"\t");
			for (int i=0;i<maps.size();i++){
				sb.append(trueHits.get(i)[j]);
				if (i==maps.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_rankedMotifOffsets_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
	}
	
	
	private void readPeakLists(){
        Vector<String> peakTags=new Vector<String>();
        for(String s : args)
        	if(s.contains("peakCall"))
        		if(!peakTags.contains(s))
        			peakTags.add(s);
    	
        // each tag represents a peak caller output file
        List<File> peakFiles=new ArrayList<File>();
        for(String tag : peakTags){
        	String name="";
        	if(tag.startsWith("--peakCall")){
        		name = tag.replaceFirst("--peakCall", ""); 
        		methodNames.add(name);
        	}
        	List<File> files = Args.parseFileHandles(args, "peakCall"+name);
        	File file= null;
        	if (!files.isEmpty()){
        		file = files.get(0);
        	}
        	peakFiles.add(file);
        }

        for (int i=0;i<methodNames.size();i++){
        	String name = methodNames.get(i);
        	String filePath = peakFiles.get(i).getAbsolutePath();
			ArrayList<Point> peakPoints = new ArrayList<Point>(); 
        	if (name.contains("GPS")){
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				// sort by descending pValue (-log P-value)
				if (!isPreSorted)
					Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
					    public int compare(GPSPeak o1, GPSPeak o2) {
					        return o1.compareByPValue(o2);
					    }
					});
				for (GPSPeak p: gpsPeaks){
					peakPoints.add(p);
				}
        	}
        	if (name.contains("SISSRS")){
        		// assume be sorted by pvalue if (!isPreSorted)
				peakPoints = load_SISSRS_File(filePath);
        	}
        	if (name.contains("MACS")){
    			List<MACSPeakRegion> macsPeaks = MACSParser.parseMACSOutput(filePath, genome);
    			if (!isPreSorted)
	    			Collections.sort(macsPeaks, new Comparator<MACSPeakRegion>(){
					    public int compare(MACSPeakRegion o1, MACSPeakRegion o2) {
					        return o1.compareByPValue(o2);
					    }
					});   			

				for (MACSPeakRegion p: macsPeaks){
    				peakPoints.add(p.getPeak());
    			}
        	}
        	if (name.contains("cisGenome")){
				peakPoints = loadCgsPointFile(filePath);
        	}  	
        	if (name.contains("QuEST")){
				peakPoints = loadCgsPointFile(filePath);
        	}  	
        	if (name.contains("PICS")){
				peakPoints = loadCgsPointFile(filePath);
        	}  	
			peaks.add(peakPoints);
        }
	}
	
	// counting number of event calls in the regions with clustered motifs
	private void proximalEventAnalysis(){
		long tic = System.currentTimeMillis();

		readPeakLists();
		
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		for (int i=0;i<peaks.size();i++){
			minCount = Math.min(minCount, peaks.get(i).size());
			maxCount = Math.max(maxCount, peaks.get(i).size());
		}
		if (rank==0){
			rank = minCount;
		}
		// print some information
		String msg = String.format("Motif score cutoff=%.2f, window size=%d, top %d peaks.\n",
				motifThreshold, windowSize, rank);
		System.out.print(msg);
		for (int i=0;i<methodNames.size();i++){
			System.out.println(methodNames.get(i)+"\t"+peaks.get(i).size());
		}
		
		// get all regions that covered by any of the peak callers, 
		// so that the motif search will be done only once
		// process by chrom to be more efficient in sorting
		ArrayList<Region> allRegions = new ArrayList<Region>();
		if (restrictRegions.isEmpty()){	// search whole genome
			for (String chrom: genome.getChromList()){
				//divide by each chrom, merge regions, then add back to the list
				ArrayList<Region> byChrom = new ArrayList<Region>();
				for (ArrayList<Point> ps:peaks){
					for (Point p: ps){
						if (p.getChrom().equalsIgnoreCase(chrom))
							byChrom.add(p.expand(windowSize));
					}
				}
				allRegions.addAll( mergeRegions(byChrom));
			}
		}
		else{	// if only search in the restrict regions
			Set<String> restrictChroms = new HashSet<String>();
			for (Region r: restrictRegions){
				restrictChroms.add(r.getChrom());
			}
			for (String chrom: genome.getChromList()){
				//divide by each chrom, merge regions, then add back to the list
				if (!restrictChroms.contains(chrom))
					continue;
				ArrayList<Region> byChrom = new ArrayList<Region>();
				for (ArrayList<Point> ps:peaks){
					for (Point p: ps){
						boolean inRestrictRegion=false;
						for (Region r: restrictRegions){
							if (r.contains(p))
								inRestrictRegion=true;
						}
						if (!inRestrictRegion)
							continue;
						if (p.getChrom().equalsIgnoreCase(chrom))
							byChrom.add(p.expand(windowSize));
					}
				}
				allRegions.addAll( mergeRegions(byChrom));
			}
		}

		System.out.println(BindingMixture.timeElapsed(tic));
				
		// get all the strong motifs in these regions
		Pair<ArrayList<Point>, ArrayList<Double>> results = getAllMotifs(allRegions, motifThreshold);
		ArrayList<Point> allMotifs = results.car();
		ArrayList<Double> allMotifScores = results.cdr();
		System.out.println(BindingMixture.timeElapsed(tic));
		
		// get all the n-ary motifs (within interDistance)
		int interDistance = 500;
		StringBuilder sb = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		for (String chrom: genome.getChromList()){
			ArrayList<ArrayList<Point>> motifClusters = new ArrayList<ArrayList<Point>>();
			
			//divide by each chrom, find cluster of motifs
			ArrayList<Point> byChrom = new ArrayList<Point>();
			for (Point motif:allMotifs){
				if (motif.getChrom().equalsIgnoreCase(chrom))
					byChrom.add(motif);
			}
			Collections.sort(byChrom);
			boolean isNewCluster = true;
			ArrayList<Point> cluster=null;
			for (int i=1; i<byChrom.size(); i++){
				if (byChrom.get(i).getLocation()-byChrom.get(i-1).getLocation()<=interDistance){
					if (isNewCluster){
						cluster = new ArrayList<Point>();
						cluster.add(byChrom.get(i-1));
						cluster.add(byChrom.get(i));
					}
					else{
						cluster.add(byChrom.get(i));
					}
					// look ahead
					if (i<byChrom.size()-1){
						if (byChrom.get(i+1).getLocation()-byChrom.get(i).getLocation()<=interDistance){
							isNewCluster = false;
						}
						else{
							motifClusters.add(cluster);
							isNewCluster = true;
						}
					}
					else{	// this is last motif
						motifClusters.add(cluster);
					}
				}
			}
			
			// print each motif clusters with event counts
			for (ArrayList<Point> tmpCluster: motifClusters){
				Region r = new Region(genome, chrom, tmpCluster.get(0).getLocation(), tmpCluster.get(tmpCluster.size()-1).getLocation());
				r=r.expand(windowSize, windowSize);
				sb.append(r.toString()).append("\t").append(tmpCluster.size());
				sb2.append(r.toString()).append("\t").append(tmpCluster.size());
				for (int i=0;i<peaks.size();i++){
					int count=0;
					ArrayList<Point> ps = peaks.get(i);
					ArrayList<Point> events = new ArrayList<Point>();
					for (int j=0;j<rank;j++){
						if (r.contains(ps.get(j))){
							events.add(ps.get(j));
							count++;
						}
					}
					sb.append("\t").append(count);
					// more stringent criteria: the events should be less than 500bp apart
					if (count>=2){
						Collections.sort(events);
						HashSet<Point> unique = new HashSet<Point>();
						for (int k=1;k<events.size();k++){
							if (events.get(k).distance(events.get(k-1))<=interDistance){
								unique.add(events.get(k));
								unique.add(events.get(k-1));
							}
						}
						count = unique.size();
					}
					sb2.append("\t").append(count);
				}
				sb.append("\n");
				sb2.append("\n");
			}
		}//each chrom
		
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_clusteredMotifEvents_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_clusteredMotifEvents500bp_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb2.toString());
	}
	
	// counting number of event calls in the regions with SINGLE motif 
	// this is to check the false positives of joint event calls
	private void singleMotifEventAnalysis(){
		long tic = System.currentTimeMillis();

		readPeakLists();
		
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		for (int i=0;i<peaks.size();i++){
			minCount = Math.min(minCount, peaks.get(i).size());
			maxCount = Math.max(maxCount, peaks.get(i).size());
		}
		if (rank==0){
			rank = minCount;
		}
		// print some information
		String msg = String.format("Motif score cutoff=%.2f, window size=%d, top %d peaks.\n",
				motifThreshold, windowSize, rank);
		System.out.print(msg);
		for (int i=0;i<methodNames.size();i++){
			System.out.println(methodNames.get(i)+"\t"+peaks.get(i).size());
		}
		
		// get all regions that covered by any of the peak callers, 
		// so that the motif search will be done only once
		// process by chrom to be more efficient in sorting
		ArrayList<Region> allRegions = new ArrayList<Region>();
		if (restrictRegions.isEmpty()){	// search whole genome
			for (String chrom: genome.getChromList()){
				//divide by each chrom, merge regions, then add back to the list
				ArrayList<Region> byChrom = new ArrayList<Region>();
				for (ArrayList<Point> ps:peaks){
					for (Point p: ps){
						if (p.getChrom().equalsIgnoreCase(chrom))
							byChrom.add(p.expand(windowSize));
					}
				}
				allRegions.addAll( mergeRegions(byChrom));
			}
		}
		else{	// if only search in the restrict regions
			Set<String> restrictChroms = new HashSet<String>();
			for (Region r: restrictRegions){
				restrictChroms.add(r.getChrom());
			}
			for (String chrom: genome.getChromList()){
				//divide by each chrom, merge regions, then add back to the list
				if (!restrictChroms.contains(chrom))
					continue;
				ArrayList<Region> byChrom = new ArrayList<Region>();
				for (ArrayList<Point> ps:peaks){
					for (Point p: ps){
						boolean inRestrictRegion=false;
						for (Region r: restrictRegions){
							if (r.contains(p))
								inRestrictRegion=true;
						}
						if (!inRestrictRegion)
							continue;
						if (p.getChrom().equalsIgnoreCase(chrom))
							byChrom.add(p.expand(windowSize));
					}
				}
				allRegions.addAll( mergeRegions(byChrom));
			}
		}

		System.out.println(BindingMixture.timeElapsed(tic));
				
		// get all the strong motifs in these regions
		Pair<ArrayList<Point>, ArrayList<Double>> results = getAllMotifs(allRegions, motifThreshold);
		ArrayList<Point> allMotifs = results.car();
		ArrayList<Double> allMotifScores = results.cdr();
		System.out.println(BindingMixture.timeElapsed(tic));
		
		// get all the single motifs (no neighbors within interDistance)
		int interDistance = 500;
		StringBuilder sb = new StringBuilder();
		StringBuilder sb2 = new StringBuilder();
		for (String chrom: genome.getChromList()){
			ArrayList<Point> singleMotifs = new ArrayList<Point>();
			
			//divide by each chrom, find cluster of motifs
			ArrayList<Point> byChrom = new ArrayList<Point>();
			for (Point motif:allMotifs){
				if (motif.getChrom().equalsIgnoreCase(chrom))
					byChrom.add(motif);
			}
			Collections.sort(byChrom);
			int size = byChrom.size();
			if (size<=2)
				continue;
			if (byChrom.get(1).getLocation()-byChrom.get(0).getLocation()>interDistance)
				singleMotifs.add(byChrom.get(0));
			for (int i=1; i<size-1; i++){
				if (byChrom.get(i).getLocation()-byChrom.get(i-1).getLocation()>interDistance &&
						byChrom.get(i+1).getLocation()-byChrom.get(i).getLocation()>interDistance	){
					singleMotifs.add(byChrom.get(i));
				}
			}
			if (byChrom.get(size-1).getLocation()-byChrom.get(size-2).getLocation()>interDistance)
				singleMotifs.add(byChrom.get(size-1));
			
			// print each single motif sites with event counts
			for (Point p: singleMotifs){
				Region r = p.expand(windowSize);
				sb.append(r.toString()).append("\t").append(0);
				sb2.append(r.toString()).append("\t").append(0);
				for (int i=0;i<peaks.size();i++){
					int count=0;
					ArrayList<Point> ps = peaks.get(i);
					ArrayList<Point> events = new ArrayList<Point>();
					for (int j=0;j<rank;j++){
						if (r.contains(ps.get(j))){
							events.add(ps.get(j));
							count++;
						}
					}
					sb.append("\t").append(count);
					// more stringent criteria: the events should be less than 500bp apart
					if (count>=2){
						Collections.sort(events);
						HashSet<Point> unique = new HashSet<Point>();
						for (int k=1;k<events.size();k++){
							if (events.get(k).distance(events.get(k-1))<=interDistance){
								unique.add(events.get(k));
								unique.add(events.get(k-1));
							}
						}
						count = unique.size();
					}
					sb2.append("\t").append(count);
				}
				sb.append("\n");
				sb2.append("\n");
			}
		}//each chrom
		
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_singleMotifEvents_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_singleMotifEvents500bp_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb2.toString());
	}
	private void getUnaryEventList(){
		long tic = System.currentTimeMillis();
		readPeakLists();
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		for (int i=0;i<peaks.size();i++){
			System.out.println(methodNames.get(i)+"\t"+peaks.get(i).size());
			minCount = Math.min(minCount, peaks.get(i).size());
			maxCount = Math.max(maxCount, peaks.get(i).size());
		}
		if (rank==0){
			rank = maxCount;		// take all possible overlaps, ranking is ignored
		}
		
		// get all regions that covered by any of the peak callers, 
		// so that the motif search will be done only once
		// process by chrom to be more efficient in sorting
		ArrayList<Region> allRegions = new ArrayList<Region>();

		for (String chrom: genome.getChromList()){
			ArrayList<Region> byChrom = new ArrayList<Region>();
			for (ArrayList<Point> ps:peaks){
				for (Point p: ps){
					if (p.getChrom().equalsIgnoreCase(chrom))
						byChrom.add(p.expand(windowSize));
				}
			}
			allRegions.addAll( mergeRegions(byChrom));
		}
		System.out.println(BindingMixture.timeElapsed(tic));
				
		// get all the strong motifs in these regions
		Pair<ArrayList<Point>, ArrayList<Double>> results = getAllMotifs(allRegions, motifThreshold);
		ArrayList<Point> allMotifs = results.car();
		ArrayList<Double> allMotifScores = results.cdr();
		System.out.println(BindingMixture.timeElapsed(tic));
		
		// Get the set of motif matches for all peak calls
		for (int i=0;i<peaks.size();i++){
			ArrayList<Point> ps = peaks.get(i);
			ArrayList<Point> toRemove = new ArrayList<Point>();
			Collections.sort(ps);
			// remove proximal events
			for (int j=1;j<ps.size();j++){
				try{
					if (ps.get(j).distance(ps.get(j-1))<=500){
						toRemove.add(ps.get(j));
						toRemove.add(ps.get(j-1));
					}
				}
				catch (IllegalArgumentException e){
					continue; // ingore events on different chrom
				}
			}
			ps.removeAll(toRemove);
			System.out.println(methodNames.get(i)+": "+ps.size());
			maps.add(getNearestMotif2(ps, allMotifs, allMotifScores));
		}
		System.out.println(BindingMixture.timeElapsed(tic));
		
		// get the intersection (motifs shared by all methods, upto top rank)
		System.out.print("\nRunning Spatial Resolution analysis ...\n");
		SetTools<Point> setTools = new SetTools<Point>();
		Set<Point> motifs_shared =  getHighRankSites(0);
		for (int i=1;i<maps.size();i++){
			Set<Point> motifs_new =  getHighRankSites(i);
			motifs_shared = setTools.intersection(motifs_shared,motifs_new);
		}
		String msg = String.format("Motif score cutoff=%.2f, window size=%d, %d shared motifs in top %d peaks.",
				motifThreshold, windowSize, motifs_shared.size(), rank);
		System.out.println(msg);
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append(msg+"\n");
		sb.append("MotifHit\tChrom\t");
		for (int i=0;i<methodNames.size();i++){
			sb.append(methodNames.get(i)+"\t");
		}
		sb.append("\n");
		for(Point motif:motifs_shared){
			sb.append(motif.toString()+"\t");
			sb.append(motif.getChrom()+"\t");
			for (int i=0;i<maps.size();i++){
				sb.append(maps.get(i).get(motif).offset);
				if (i==maps.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		BindingMixture.writeFile(outName+"_"+methodNames.size()+"methods_unaryEventLists_"
				+String.format("%.2f_",motifThreshold)
				+windowSize+".txt", sb.toString());
		
	}
	// merge the overlapped regions 
	// if "toExpandRegion"=true, expand each region on both side to leave enough space, 
	// to include every potential reads, then merge
	private ArrayList<Region> mergeRegions(ArrayList<Region> regions){
		ArrayList<Region> mergedRegions = new ArrayList<Region>();
		if (regions.isEmpty())
			return mergedRegions;
		Collections.sort(regions);
		Region previous = regions.get(0);
		mergedRegions.add(previous);
		
		for (Region region: regions){
			// if overlaps with previous region, combine the regions
			if (previous.overlaps(region)){
				mergedRegions.remove(previous);
				previous = previous.combine(region);
			}
			else{
				previous = region;
			}
			mergedRegions.add(previous);
		}
		return mergedRegions;
	}//end of mergeRegions method
	
	private Set<Point> getHighRankSites(int method){
		HashMap<Point, MotifHit> m = maps.get(method);
		Set<Point> motifs = m.keySet();
		Set<Point> motifs_top = new HashSet<Point>();
		maps.get(0).keySet();
		for (Point p: motifs){
			if (m.get(p).rank < rank)
				motifs_top.add(p);
		}
		return motifs_top;
	}

	// get the all nearest motif hits in the region (windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with the nearest peak
	private HashMap<Point, MotifHit> getNearestMotif(ArrayList<Point> peaks, int rank, double threshold){	
		HashMap<Point, MotifHit> motifs = new HashMap<Point, MotifHit>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		int range = Math.min(peaks.size(), rank); // only compare top #(rank) peaks
		for(int i=0; i<range; i++){
			if (i % 1000==0)
				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(windowSize);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			//search for whole region
			for(int z=0; z<r.getWidth(); z++){		
				double score = profiler.getMaxScore(z);
				if(score >= threshold){
					Point motifPos = new Point(genome, peak.getChrom(), r.getStart()+z+profiler.getMatrix().length()/2);
					MotifHit hit = new MotifHit(motifPos, score, peak, i);
					//if other peaks also associated with this motif, keep the nearest one
					if (motifs.containsKey(motifPos)){
						if (Math.abs(hit.offset)<Math.abs(motifs.get(motifPos).offset))
							motifs.put(motifPos, hit);
					}
					else{
						motifs.put(motifPos, hit);
					}
				}
			}
		}
		return motifs;
	}
	// get the all nearest motif hits in the region (windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with the nearest peak
	private HashMap<Point, MotifHit> getNearestMotif2(ArrayList<Point> peaks, ArrayList<Point> allMotifs, ArrayList<Double>allMotifScores){	
		HashMap<Point, MotifHit> motifs = new HashMap<Point, MotifHit>();
		int numPeaksWithinMotifs = 0;
		for(int i=0;i<peaks.size();i++){
			Point peak =peaks.get(i);
			int nearestIndex = -1;
			int nearestDistance = windowSize;
			//TODO: should search for the position, not iterate everyone
			for(int j=0;j<allMotifs.size();j++){
				Point motif =allMotifs.get(j);
				if (peak.getChrom().equalsIgnoreCase(motif.getChrom())){
					int distance = peak.distance(motif);
					if (distance<=nearestDistance){
						nearestDistance = distance;
						nearestIndex = j;
					}
				}
			}
			
			if (nearestIndex !=-1){			// motif hit is within the window
				numPeaksWithinMotifs++;
				Point nearestMotif = allMotifs.get(nearestIndex);
				MotifHit hit = new MotifHit(nearestMotif, allMotifScores.get(nearestIndex), peak, i);
				//if other peaks also associated with this motif, keep the nearest one
				if (motifs.containsKey(nearestMotif)){
					if (Math.abs(hit.offset)<Math.abs(motifs.get(nearestMotif).offset))
						motifs.put(nearestMotif, hit);
				}
				else{
					motifs.put(nearestMotif, hit);
				}
			}
		}
		System.out.printf("%d/%d (%.1f%%)%n", numPeaksWithinMotifs, peaks.size(), 100.0*((double)numPeaksWithinMotifs)/(double)peaks.size(), windowSize);
		return motifs;
	}
	
	// this is an alternative test to exclude the possibility that one method may 
	// take advantage by predicting multiple peaks around the motif position
	// therefore, if a motif have more than one peaks nearby, do not include it
	private HashMap<Point, MotifHit> getNearestMotif3(ArrayList<Point> peaks, ArrayList<Point> allMotifs, ArrayList<Double>allMotifScores){	
		HashMap<Point, MotifHit> motifs = new HashMap<Point, MotifHit>();
		for(int j=0;j<allMotifs.size();j++){
			Point motif =allMotifs.get(j);
			int count = 0;
			int peakIndex = -1;
			for(int i=0;i<peaks.size();i++){
				Point peak =peaks.get(i);
				int nearestDistance = windowSize;
				if (peak.getChrom().equalsIgnoreCase(motif.getChrom())){
					int distance = peak.distance(motif);
					if (distance<=nearestDistance){
						count++;
						nearestDistance = distance;
						peakIndex = i;
					}
				}
			}
			if (count!=1)	// not exactly one peak in the 100bp to the motif, exclude this motif
				continue;
			
			Point nearestMotif = allMotifs.get(j);
			MotifHit hit = new MotifHit(nearestMotif, allMotifScores.get(j), peaks.get(peakIndex), peakIndex);
			motifs.put(nearestMotif, hit);
		}
		return motifs;
	}
		
	// get the all nearest motif hits in the region (windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with the nearest peak
	private Pair<ArrayList<Point>, ArrayList<Double>> getAllMotifs(ArrayList<Region> regions, double threshold){
		ArrayList<Point> allMotifs = new ArrayList<Point>();
		ArrayList<Double> scores = new ArrayList<Double>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		int count = regions.size();
		System.out.print("\nGetting motifs from "+count+" regions ...\n");
		for(int i=0; i<count; i++){
			if (i % 1000==0)
				System.out.println(i);
			Region r= regions.get(i);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			//search for whole region
			for(int z=0; z<r.getWidth(); z++){		
				double score = profiler.getMaxScore(z);
				if(score >= threshold){
					Point motifPos = new Point(genome, r.getChrom(), r.getStart()+z+profiler.getMatrix().length()/2);
					allMotifs.add(motifPos);
					scores.add(score * (profiler.getMaxStrand(z)=='+'?1:-1)); // positive score if motif match on '+' strand
				}
			}
		}
		return new Pair<ArrayList<Point>, ArrayList<Double>>(allMotifs, scores);
	}
	
	private class MotifHit{
		Point motif;
		double score;
		Point peak;
		int rank;	// rank of the peak in the prediction list, method dependent
		int offset;	// positive: peak is on 5' direction of motif, negative: otherwise
		MotifHit(Point motif, double score, Point peak, int rank){
			this.motif = motif;
			this.score = score;
			this.peak = peak;
			this.rank = rank;
			offset = peak.offset(motif)*(score>0?1:-1);
		}
	}

	// load text file in CGS Point format
	// chr:coord, e.g. 1:234234
	private ArrayList<Point> loadCgsPointFile(String filename) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				Region point = Region.fromString(genome, line);
				points.add(new Point(genome, point.getChrom(),point.getStart()));
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
		return points;
	}

	/** load text file in SISSRS output BED format, then sort by p-value
	 * 	Chr		cStart	cEnd	NumTags	Fold	p-value
		---		------	----	-------	----	-------
		chr1	4132791	4132851	38		71.44	5.0e-006
	 */
	private ArrayList<Point> load_SISSRS_File(String filename) {

		File file = new File(filename);
		FileReader in = null;
		BufferedReader bin = null;
		ArrayList<Point> points = new ArrayList<Point>();
		ArrayList<SISSRS_Peak> peaks = new ArrayList<SISSRS_Peak>();
		try {
			in = new FileReader(file);
			bin = new BufferedReader(in);
			// skip header, data starts at 58th line
			for (int i=1;i<58;i++){
				bin.readLine();
			}
			String line;
			while((line = bin.readLine()) != null) { 
				line = line.trim();
				String[] t = line.split("\\t");
				if (t.length==6){	// region format
					peaks.add(new SISSRS_Peak(genome, t[0].replaceFirst("chr", ""), Integer.parseInt(t[1]), Integer.parseInt(t[2]),
						Double.parseDouble(t[3]), Double.parseDouble(t[4]), Double.parseDouble(t[5])));
				}
				if (t.length==5){	// point format
					peaks.add(new SISSRS_Peak(genome, t[0].replaceFirst("chr", ""), Integer.parseInt(t[1]), Integer.parseInt(t[1]),
						Double.parseDouble(t[2]), Double.parseDouble(t[3]), Double.parseDouble(t[4])));
				}
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
		// sort by pvalue
		if (!isPreSorted)
			Collections.sort(peaks);
		
		for (SISSRS_Peak p:peaks){
			points.add(p.getPeak());
		}
		return points;
	}
	///////////////////////////////////////////////////////////////
	// old code
	
	private void printMotifOffsetDistance(){
//		if (peaks_macs.size()>0)
//			motifOffsetDistance(peaks_macs, "MACS");
		if (peaks_gps.size()>0)
			motifOffsetDistance(peaks_gps, "GPS");
	}
	
	// here we only print out the offset distance between motif and peak
	// we use another matlab code to calculate 
	// 1. the percentage of motif occurence
	// 2. the spatial resolution
	private void motifOffsetDistance(ArrayList<Point> peaks, String method){
		System.out.print("\nRunning Motif Occurrence analysis ...\n");
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);

		// for each peak and motif threshold, the shortest distance between BS and motif
		int[][]distance = new int[motifThresholds.length][peaks.size()];
		for(int i=0; i<peaks.size(); i++){
			if (i % 100==0)
				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(windowSize);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			for(int j=0;j<motifThresholds.length;j++){
				double threshold = motifThresholds[j];
				//search from BS outwards
				for(int z=0; z<=r.getWidth()/2; z++){
					double leftScore= profiler.getMaxScore(windowSize-z);
					double rightScore= profiler.getMaxScore(windowSize+z);				
					if(rightScore>=threshold){
						distance[j][i] = z;
						break;
					}
					if(leftScore>=threshold){
						distance[j][i] = -z;
						break;
					}
					// if motif score at this position is too small, and reach end of region
					if (z==r.getWidth()/2){
						distance[j][i] = NOHIT_OFFSET;
					}
				}
			}
		}
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+" "+method+"\n");
		sb.append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++){
			sb.append("\t"+motifThresholds[j]);
		}
		sb.append("\n");
		for(int i=0; i<peaks.size(); i++){
			sb.append(peaks.get(i).getLocationString());
			for(int j=0;j<motifThresholds.length;j++){
				sb.append("\t"+distance[j][i]);
			}
			sb.append("\n");
		}
		BindingMixture.writeFile("CTCF_motif_distance_byScore_"+windowSize+"_"+method+".txt", sb.toString());
	}	
		
	private void motifBasedAnalysis(){
		long tic = System.currentTimeMillis();
		// Get the set of motif matches from all methods
		System.out.print("\nGetting motifs from GPS result ...\n");
		maps_gps = getAllNearestMotifs(peaks_gps, motifThresholds[0]);
		System.out.print("\nGetting motifs from MACS result ...\n");
		maps_macs = getAllNearestMotifs(peaks_macs, motifThresholds[0]);
		System.out.println(BindingMixture.timeElapsed(tic));
		sharedMotif_SpatialResolution();
//		bindingSpatialResolutionHistrogram();
//		motifCoverage();
	}
	
	// for each shared motif, average all the peaks that fits the criteria
	private void sharedMotif_SpatialResolution(){
		System.out.print("\nRunning Spatial Resolution analysis ...\n");
		long tic = System.currentTimeMillis();
		SetTools<Point> setTools = new SetTools<Point>();
		// for each peak and motif threshold, the shortest distance between BS and motif
		int peakCount = Math.min(peaks_gps.size(), peaks_macs.size());
		int step = 100;
		int stepCount = peakCount/step;

		// array to store average spatial resolution, [diff methods][motif scores][peak rank steps]
		double[][][]resolution = new double[2][motifThresholds.length][stepCount];

		for(int i=0;i<motifThresholds.length;i++){
			// get only the motifs that pass the score cutoff
			Set<Point> motifs_gps = new HashSet<Point>();
			Set<Point> motifs_macs = new HashSet<Point>();
			for (Point motif: maps_gps.keySet())
				if (maps_gps.get(motif).get(0).score >= motifThresholds[i])
					motifs_gps.add(motif);
			for (Point motif: maps_macs.keySet())
				if (maps_macs.get(motif).get(0).score >= motifThresholds[i])
					motifs_macs.add(motif);
			
			// intersection set for spatial resolution comparison
			// only consider the motif positions that are cover by all methods
			Set<Point> intersection = setTools.intersection(motifs_gps,motifs_macs);
			System.out.print(String.format("Motif score cutoff=%.2f,\t%d\t shared motifs\n",
					motifThresholds[i], intersection.size()));
			
			// calculate the average spatial resolution (offset) in the range of rank steps
			for (int j=0;j<stepCount;j++) {
				int total = 0;
				int count = 0;
				int total_macs = 0;
				int count_macs = 0;
				for (Point motif:intersection){				// for each shared motif
					for (MotifHit hit:maps_gps.get(motif))
						if (hit.rank<=(j+1)*step){			// if the peak is in top rank range
							total+= Math.abs(hit.offset);
							count++;
						}
					for (MotifHit hit:maps_macs.get(motif))
						if (hit.rank<=(j+1)*step){
							total_macs+= Math.abs(hit.offset);
							count_macs++;
						}
				}
				resolution[0][i][j] = total/(double)count;
				resolution[1][i][j] = total_macs/(double)count_macs;
			}
		}
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n");
		sb.append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%.2f", motifThresholds[j]));
		sb.append("\n");
		
		for(int j=0; j<stepCount; j++){
			sb.append((j+1)*step);
			for(int i=0;i<motifThresholds.length;i++)
				sb.append("\t"+resolution[0][i][j]);
			sb.append("\n");
		}
		BindingMixture.writeFile("CTCF_GPS_SpatialResolution_"+windowSize+".txt", sb.toString());

		sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n");
		sb.append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%.2f", motifThresholds[j]));
		sb.append("\n");
		
		for(int j=0; j<stepCount; j++){
			sb.append((j+1)*step);
			for(int i=0;i<motifThresholds.length;i++)
				sb.append("\t"+resolution[1][i][j]);
			sb.append("\n");
		}
		BindingMixture.writeFile("CTCF_MACS_SpatialResolution_"+windowSize+".txt", sb.toString());
		System.out.println(BindingMixture.timeElapsed(tic));
	}
	
	// for all the shared motif, pick the nearest peaks
	// for all kinds of threshold
	private void bindingSpatialResolutionHistogram(){
		System.out.print("\nRunning Spatial Resolution analysis ...\n");
		SetTools<Point> setTools = new SetTools<Point>();
		// for each peak and motif threshold, the shortest distance between BS and motif
		for(int m=0;m<motifThresholds.length;m++){
			// get only the motifs that pass the score cutoff
			Set<Point> motifs_gps = new HashSet<Point>();
			Set<Point> motifs_macs = new HashSet<Point>();
			for (Point motif: maps_gps.keySet())
				if (maps_gps.get(motif).get(0).score >= motifThresholds[m])
					motifs_gps.add(motif);
			for (Point motif: maps_macs.keySet())
				if (maps_macs.get(motif).get(0).score >= motifThresholds[m])
					motifs_macs.add(motif);
			
			// intersection set for spatial resolution comparison
			// only consider the motif positions that are cover by all methods
			ArrayList<Point> sharedMotifs = new ArrayList<Point>();
			sharedMotifs.addAll( setTools.intersection(motifs_gps,motifs_macs));
			System.out.print(String.format("Motif score cutoff=%.2f,\t%d\t shared motifs\n",
					motifThresholds[m], sharedMotifs.size()));
			// calculate the spatial resolution (offset) 
			int rank = 5000;
			// array to store all spatial resolution values, [methods][shared peaks with motif]
			double[][]resolution = new double[2][sharedMotifs.size()];
			for (int i=0;i<sharedMotifs.size();i++){	// for each shared motif
				Point motif = sharedMotifs.get(i);
				int minDistance = Integer.MAX_VALUE;
				int minOffset=Integer.MAX_VALUE;
				for (MotifHit hit:maps_gps.get(motif)){
					if (hit.rank<=rank){			// if the peak is in top rank range
						if (minDistance > Math.abs(hit.offset)){
							minDistance = Math.abs(hit.offset);
							minOffset = hit.offset;
						}
					}
				}
				resolution[0][i]=minOffset;
				
				minDistance = Integer.MAX_VALUE;
				minOffset=Integer.MAX_VALUE;
				for (MotifHit hit:maps_macs.get(motif)){
					if (hit.rank<=rank){			// if the peak is in top rank range
						if (minDistance > Math.abs(hit.offset)){
							minDistance = Math.abs(hit.offset);
							minOffset = hit.offset;
						}
					}
				}
				resolution[1][i]=minOffset;
			}
			
			// output results
			StringBuilder sb = new StringBuilder();
			for(int j=0;j<sharedMotifs.size();j++){
				if (resolution[0][j]<=windowSize && resolution[1][j]<=windowSize){
					sb.append(sharedMotifs.get(j).toString()+"\t");
					sb.append(sharedMotifs.get(j).getChrom()+"\t");
					sb.append(String.format("%.2f", resolution[0][j])+"\t"+
							String.format("%.2f", resolution[1][j])+"\n");
				}
			}
			BindingMixture.writeFile("CTCF_SpatialResolutionHistogram_"
					+String.format("%.2f_",motifThresholds[m])
					+windowSize+".txt", sb.toString());
		}
	}
	

	// Motif Coverage analysis is to assess the sensitivity of the peak calling methods
	// We generate a high confident "true" binding motif list by all the strong motifs that
	// are covered by any of two methods, or covered by majority of 3+ methods.
	// Then for each motif score threshold, we determine how many of the motifs are covered 
	// by the top ranking peaks called by each method
	private void motifCoverage(){
		System.out.print("\nRunning Motif Coverage analysis ...\n");
		long tic = System.currentTimeMillis();
		SetTools<Point> setTools = new SetTools<Point>();
		
		int peakCount = Math.min(peaks_gps.size(), peaks_macs.size());
		int step = 100;
		int stepCount = peakCount/step;

		// array to store motif coverage count, [diff methods][motif scores][peak rank steps]
		int[][][]motifCovered = new int[2][motifThresholds.length][stepCount];
		int motifTotal[] = new int[motifThresholds.length];
		
		for(int i=0;i<motifThresholds.length;i++){
			// get only the motifs that pass the score cutoff
			Set<Point> motifs_gps = new HashSet<Point>();
			Set<Point> motifs_macs = new HashSet<Point>();
			for (Point motif: maps_gps.keySet())
				if (maps_gps.get(motif).get(0).score >= motifThresholds[i])
					motifs_gps.add(motif);
			for (Point motif: maps_macs.keySet())
				if (maps_macs.get(motif).get(0).score >= motifThresholds[i])
					motifs_macs.add(motif);
					
			//TODO
			// union set for sensitivity analysis
			// all the possible motif positions that are covered by any methods
			// if we have more methods, covered by majority of methods ...
			Set<Point> union = setTools.union(motifs_gps,motifs_macs);
			motifTotal[i]=union.size();
			System.out.print(String.format("Motif score cutoff=%.2f,\t%d\t total motifs\n",
					motifThresholds[i], union.size()));
			// calculate the average motif coverage in the range of rank steps
			for (int j=0;j<stepCount;j++) {
				int count = 0;
				int count_macs = 0;
				for (Point motif:union){				// for every motif
					if (maps_gps.containsKey(motif))
						for (MotifHit hit:maps_gps.get(motif))
							if (hit.rank<=(j+1)*step){			// if the peak is in top rank range
								count++;
								break;							// only one is enough
							}
					if (maps_macs.containsKey(motif))
						for (MotifHit hit:maps_macs.get(motif))
							if (hit.rank<=(j+1)*step){
								count_macs++;
								break;							
							}
				}
				motifCovered[0][i][j] = count;
				motifCovered[1][i][j] = count_macs;
			}
		}
		
		// output results
		StringBuilder sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n");
		sb.append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%.2f", motifThresholds[j]));
		sb.append("\n");		
		sb.append("#Total Motif");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%d", motifTotal[j]));
		sb.append("\n");
		
		for(int j=0; j<stepCount; j++){
			sb.append((j+1)*step);
			for(int i=0;i<motifThresholds.length;i++)
				sb.append("\t"+motifCovered[0][i][j]);
			sb.append("\n");
		}
		BindingMixture.writeFile("CTCF_GPS_MotifCoverage_"+windowSize+".txt", sb.toString());

		sb = new StringBuilder();
		sb.append("#"+motif.name+" "+ motif.version+"\n");
		sb.append("#Motif Score");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%.2f", motifThresholds[j]));
		sb.append("\n");		
		sb.append("#Total Motif");
		for(int j=0;j<motifThresholds.length;j++)
			sb.append("\t"+String.format("%d", motifTotal[j]));
		sb.append("\n");
		
		for(int j=0; j<stepCount; j++){
			sb.append((j+1)*step);
			for(int i=0;i<motifThresholds.length;i++)
				sb.append("\t"+motifCovered[1][i][j]);
			sb.append("\n");
		}
		BindingMixture.writeFile("CTCF_MACS_MotifCoverage_"+windowSize+".txt", sb.toString());
		System.out.println(BindingMixture.timeElapsed(tic));
	}	
	
	// get the all nearest motif hits in the region (within windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with all the peaks
	
	// returns a map of motif -- map to all associated peaks
	private HashMap<Point, ArrayList<MotifHit>> getAllNearestMotifs(ArrayList<Point> peaks, double threshold){	
		HashMap<Point, ArrayList<MotifHit>> motifs = new HashMap<Point, ArrayList<MotifHit>>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		
		for(int i=0; i<peaks.size(); i++){
			if (i % 1000==0)
				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(windowSize);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			//search for whole region
			for(int z=0; z<r.getWidth(); z++){		
				double score = profiler.getMaxScore(z);
				if(score >= threshold){
					Point motifPos = new Point(genome, peak.getChrom(), r.getStart()+z+profiler.getMatrix().length()/2);
					MotifHit hit = new MotifHit(motifPos, score, peak, i);
					//if other peaks also associate with this motif position, add to the list
					if (motifs.containsKey(motifPos))	
						motifs.get(motifPos).add(hit) ;
					else{
						ArrayList<MotifHit> hits = new ArrayList<MotifHit>();
						hits.add(hit);
						motifs.put(motifPos, hits);
					}
				}
			}
		}
		return motifs;
	}
	
	class SISSRS_Peak implements Comparable<SISSRS_Peak>{
		Region region;
		double tags;
		double fold;
		double pvalue;
		SISSRS_Peak(Genome g, String chrom, int start, int end, double tags, double fold, double pvalue){
			this.region = new Region(g, chrom, start, end);
			this.tags = tags;
			this.fold = fold;
			this.pvalue = pvalue;
		}
		Point getPeak(){
			return region.getMidpoint();
		}
		//Comparable default method
		public int compareTo(SISSRS_Peak p) {
			double diff = pvalue-p.pvalue;
			return diff==0?0:(diff<0)?-1:1;
		}
	}
}