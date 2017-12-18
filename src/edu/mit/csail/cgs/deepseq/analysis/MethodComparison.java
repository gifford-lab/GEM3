package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
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
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils.SISSRS_Event;
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

public class MethodComparison {
	private final int NOHIT_OFFSET = 999;
	
	private boolean isPreSorted = false;	
	private int windowSize = 50;	
	private int rank = 0;
	private boolean sortByStrength = false;
	
	private Genome genome;
	private String[] args;
	private WeightMatrix motif = null;
	private String outName="out";
	
	private String truePointPath = null;
	// each element in the list is for one ChIP-Seq method
	private ArrayList<String> methodNames = new ArrayList<String>();
	private ArrayList<ArrayList<Point>> events = new ArrayList<ArrayList<Point>>();
	private ArrayList<HashMap<Point, TrueHit>> maps = new ArrayList<HashMap<Point, TrueHit>>();	
	
	private ArrayList<Point> peaks_gps = new ArrayList<Point>();
	private ArrayList<Point> peaks_macs = new ArrayList<Point>();
	HashMap<Point, ArrayList<TrueHit>> maps_gps = null;
	HashMap<Point, ArrayList<TrueHit>> maps_macs = null;
	HashMap<Point, TrueHit> map_gps = null;
	HashMap<Point, TrueHit> map_macs = null;
	
	// motif score cutoffs
	//	double[] motifThresholds = new double[]{0,4,8,12,16,20,24,28};
	double[] motifThresholds = new double[]{5.587944030761719, 7.559912204742432, 
		11.52034854888916, 12.886543273925781,	15.858697891235352}; //Ctcf mouse
	private double motifThreshold=0;	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		MethodComparison analysis = new MethodComparison(args);
		// peak based analysis
		// We only print out the ChIP-exo point (at diff threshold) distances for each peak
		// spatialResolution or motifOccurence are analysis in Matlab 
		int type = Args.parseInteger(args, "analysisType", 0);
		switch(type){
			case 0:analysis.printPointOffsets();break;
			default: System.err.println("Unrecognize analysis type: "+type);
		}
	}
	
	MethodComparison(String[] args){
		this.args = args;
		ArgParser ap = new ArgParser(args);
		Set<String> flags = Args.parseFlags(args);		
		genome = CommonUtils.parseGenome(args);
	    	    
		// some parameters
		windowSize = Args.parseInteger(args, "windowSize", 50);
		rank = Args.parseInteger(args, "rank", 0);
		outName = Args.parseString(args,"out",outName);
		truePointPath = Args.parseString(args,"truePoints",truePointPath);
		sortByStrength = flags.contains("ss");	// only for GPS
		isPreSorted = flags.contains("sorted");
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
	private void printPointOffsets() throws IOException {
		long tic = System.currentTimeMillis();
		readPeakLists();
		
		int minCount = Integer.MAX_VALUE;
		int maxCount = 0;
		for (int i=0;i<methodNames.size();i++){
			minCount = Math.min(minCount, events.get(i).size());
			maxCount = Math.max(maxCount, events.get(i).size());
		}
		if (rank==0){
			rank = minCount;
		}
		// print some information
		String msg = String.format("Window size=%d, top %d events.\n",
				windowSize, rank);
		System.out.print(msg);
		System.out.println("\nNumber of events:");
		for (int i=0;i<methodNames.size();i++){
			System.out.println(methodNames.get(i)+"\t"+events.get(i).size());
		}
		
		// get all regions that covered by any of the peak callers, 
		ArrayList<Region> allRegions = new ArrayList<Region>();
		for (String chrom: genome.getChromList()){
			ArrayList<Region> byChrom = new ArrayList<Region>();
			for (ArrayList<Point> ps:events){
				for (Point p: ps){
					if (p.getChrom().equalsIgnoreCase(chrom))
						byChrom.add(p.expand(windowSize));
				}
			}
			allRegions.addAll( mergeRegions(byChrom));
		}
		System.out.println(CommonUtils.timeElapsed(tic));
				
		// get all the hits in these regions
		ArrayList<Point> allPoints = readTruePoints(truePointPath);
		
		SetTools<Point> setTools = new SetTools<Point>();
		
		System.out.printf("%n%d true hits (in %d regions).%n%n", allPoints.size(), allRegions.size());

		// Get the set of hits for all peak calls		
		System.out.printf("Events within a %d bp window to a hit:%n", windowSize);
		for (int i=0;i<events.size();i++){
			System.out.printf("%s \t#events: ", methodNames.get(i));
			HashMap<Point, TrueHit> event2trueHits = getNearestPoint(events.get(i), allPoints);
			maps.add(event2trueHits);
			System.out.printf("\t#trueHits: %d", event2trueHits.keySet().size());
			System.out.println();
		}
		System.out.println(CommonUtils.timeElapsed(tic));
		
		System.out.printf("\nHits covered by an event:\n");
		for(int i = 0; i < methodNames.size(); i++)
			System.out.printf("%s\t%d/%d (%.1f%%)\n", methodNames.get(i), maps.get(i).size(), allPoints.size(), 100.0*((double)maps.get(i).size())/(double)allPoints.size());
				
		// get the intersection (hits shared by all methods, upto top rank)
		System.out.print("\nRunning Spatial Resolution analysis ...\n");
		Set<Point> hits_shared =  getHighRankSites(0);
		for (int i=1;i<maps.size();i++){
			Set<Point> hits_new =  getHighRankSites(i);
			hits_shared = setTools.intersection(hits_shared,hits_new);
		}
		msg = String.format("Window size=%d, %d shared hits in top %d events.",
				windowSize, hits_shared.size(), rank);
		System.out.println(msg);
		
		StringBuilder args_sb = new StringBuilder();
		for (String arg:args){
			args_sb.append(arg).append(" ");
		}
		String args_str = args_sb.toString();
		
		// output results, the spatial resolution (offset) 
		StringBuilder sb = new StringBuilder();
		sb.append(args_str+"\t"+msg+"\n");
		sb.append("TrueHit\tChrom");
		for (int i=0;i<methodNames.size();i++){
			sb.append("\t"+methodNames.get(i));
		}
		sb.append("\n");
		for(Point hit:hits_shared){
			sb.append(hit.toString()+"\t");
			sb.append(hit.getChrom()+"\t");
			for (int i=0;i<maps.size();i++){
				sb.append(maps.get(i).get(hit).offset);
				if (i==maps.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		CommonUtils.writeFile(outName+"_"+methodNames.size()+"methods_sharedPointOffsets_"
				+windowSize+".txt", sb.toString());
		
// ************************************************************************************
// ************************************************************************************
		
		// get the union (hits at least by one of all methods)
		Set<Point> hits_union = maps.get(0).keySet();
		for (int i=1;i<maps.size();i++){
			Set<Point> hits_new = maps.get(i).keySet();
			hits_union = setTools.union(hits_union,hits_new);
		}
		msg = String.format("Total %d hits from all methods.", hits_union.size());
		System.out.println(msg);
		
		// output results, the spatial resolution (offset) 
		sb = new StringBuilder();
		sb.append(args_str+"\t"+msg+"\n");
		sb.append("TrueHit\tChrom\t");
		for (int i=0;i<methodNames.size();i++){
			sb.append(methodNames.get(i)+"_offset\t");
			sb.append(methodNames.get(i)+"_rank\t");
		}
		sb.append("\n");
		for(Point hit:hits_union){
			sb.append(hit.toString()+"\t");
			sb.append(hit.getChrom()+"\t");
			for (int i=0;i<maps.size();i++){
				HashMap<Point, TrueHit> m = maps.get(i);
				if (m.containsKey(hit)){
					sb.append(m.get(hit).offset).append("\t");
					sb.append(m.get(hit).rank);
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
		CommonUtils.writeFile(outName+"_"+methodNames.size()+"methods_allPointOffsets_"
				+windowSize+".txt", sb.toString());
		
		System.out.println();
		
// ************************************************************************************
// ************************************************************************************
		// motif occurrence / motif coverage (sensitivity)
		// get the motif offset of the peaks (rank indexed by each method)
		// if no motif hit, offset=NOHIT=999
		// print out the offset, for Matlab processing and plotting
		
		System.out.println("Get ranked peak-hit offset list ... ");
		
		ArrayList<HashMap<Point, Integer>> allPeakOffsets = new ArrayList<HashMap<Point, Integer>>();
		for (int i=0; i<methodNames.size();i++){
			allPeakOffsets.add(peak2HitOffset(events.get(i), allPoints));
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
			for (int i=0;i<allPeakOffsets.size();i++){
				if (j<events.get(i).size()){
					Point peak = events.get(i).get(j);
					sb.append(allPeakOffsets.get(i).containsKey(peak)?allPeakOffsets.get(i).get(peak):NOHIT_OFFSET);
				}
				else
					sb.append(NOHIT_OFFSET);
				
				if (i==allPeakOffsets.size()-1)
					sb.append("\n");
				else
					sb.append("\t");
			}
		}
		CommonUtils.writeFile(outName+"_"+methodNames.size()+"methods_rankedPointOffsets_"
				+windowSize+".txt", sb.toString());
		System.out.println("Done! " + CommonUtils.timeElapsed(tic));
	}

	private void readPeakLists() throws IOException {
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
        	if (name.contains("GPS")||name.contains("GEM")){
				List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
				// sort by descending pValue (-log P-value)
				if (!isPreSorted){
					if (sortByStrength){
						Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
						    public int compare(GPSPeak o1, GPSPeak o2) {
						        return o1.compareByIPStrength(o2);
						    }
						});	
					}
					else{
						Collections.sort(gpsPeaks, new Comparator<GPSPeak>(){
						    public int compare(GPSPeak o1, GPSPeak o2) {
						        return o1.compareByPV_lg10(o2);
						    }
						});
					}
				}
				for (GPSPeak p: gpsPeaks){
					peakPoints.add(p);
				}
        	}
        	else if (name.contains("SISSRS")){
        		// assume be sorted by pvalue if (!isPreSorted)
        		ArrayList<SISSRS_Event> SISSRs_events = CommonUtils.load_SISSRs_events(genome, filePath, isPreSorted);
        		for (SISSRS_Event se: SISSRs_events){
        			peakPoints.add(se.getPeak());
        		}
        	}
        	else if (name.contains("MACS")){
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
        	else{
				peakPoints = CommonUtils.loadCgsPointFile(filePath, genome);
        	}  

        	events.add(peakPoints);
        }
	}
	
	// read true peak points and sort by location
	private ArrayList<Point> readTruePoints (String filePath) throws IOException {
		ArrayList<Point> truePoints = new ArrayList<Point>();
		List<GPSPeak> gpsPeaks = GPSParser.parseGPSOutput(filePath, genome);
		// sort by location
		Collections.sort(gpsPeaks);
		for (GPSPeak p: gpsPeaks){
			truePoints.add(p);
		}
		return truePoints;
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
	
	// get the top ranking sites of a set of calls
	private Set<Point> getHighRankSites(int method){
		HashMap<Point, TrueHit> m = maps.get(method);
		Set<Point> motifs = m.keySet();
		Set<Point> motifs_top = new HashSet<Point>();
		maps.get(0).keySet();
		for (Point p: motifs){
			if (m.get(p).rank < rank)
				motifs_top.add(p);
		}
		return motifs_top;
	}
	
	// get the nearest true peakPoint in the region (windowSize) around each event in the list
	// this is peakPoint-centered, some event may associate with a peakPoint but then took away by other peaks
	// if a peakPoint have more than one events nearby, associate it with the nearest event
	// Assuming the peakPoint positions are sorted
	private HashMap<Point, TrueHit> getNearestPoint(ArrayList<Point> events, ArrayList<Point> allPoints){
		// make a copy and sort
		ArrayList<Point> ps = (ArrayList<Point>) events.clone();
		Collections.sort(ps);
		
		HashMap<Point, TrueHit> trueHits = new HashMap<Point, TrueHit>();
		int numPeaksWithinMotifs = 0;
		int firstMatchIndex = 0;		// the index of first motif match index for previous peak
										// this will be use as start search position of the inner loop
		for(int i=0;i<ps.size();i++){
			Point peak =ps.get(i);
			int nearestIndex = -1;
			int nearestDistance = windowSize;
			boolean firstMatchFound = false;
			// search for the position, not iterate everyone
			for(int j=firstMatchIndex;j<allPoints.size();j++){
				Point motif =allPoints.get(j);
				if (peak.getChrom().equalsIgnoreCase(motif.getChrom())){
					int distance = peak.distance(motif);
					if (distance<=nearestDistance){
						if (!firstMatchFound){		// update on the first match
							firstMatchIndex = j;
							firstMatchFound = true;
						}							
						nearestDistance = distance;
						nearestIndex = j;
					}
					else{		// if distance > nearest, and found a match already, it is getting farther away, stop search 
						if (firstMatchFound){	
							break;
						}
					}					
				}
			}
			
			if (nearestIndex !=-1){			// hit is within the window
				numPeaksWithinMotifs++;
				Point nearestPoint = allPoints.get(nearestIndex);
				TrueHit hit = new TrueHit(nearestPoint, peak, i);
				//if other peaks also associated with this motif, keep the nearest one
				if (trueHits.containsKey(nearestPoint)){
					if (Math.abs(hit.offset)<Math.abs(trueHits.get(nearestPoint).offset))
						trueHits.put(nearestPoint, hit);
				}
				else{
					trueHits.put(nearestPoint, hit);
				}
			}
		}
		System.out.printf("%d/%d (%.1f%%)", numPeaksWithinMotifs, events.size(), 100.0*((double)numPeaksWithinMotifs)/(double)events.size(), windowSize);
		return trueHits;
	}
	
	// get the nearest motif offset in the region (windowSize) around each peak in the list
	// this is peak-centered, all peaks will have a offset, NOHIT_OFFSET for no motif found.
	private HashMap<Point, Integer> peak2HitOffset(ArrayList<Point> peaks, ArrayList<Point> allMotifs){
		// make a copy of the list and sort
		ArrayList<Point> ps = (ArrayList<Point>) peaks.clone();
		Collections.sort(ps);
		
		HashMap<Point, Integer> peaksOffsets = new HashMap<Point, Integer>();
		int firstMatchIndex = 0;		// the index of first motif match index for previous peak
										// this will be use as start search position of the inner loop
		for(int i=0;i<ps.size();i++){
			Point peak =ps.get(i);
			int nearestIndex = -1;
			int nearestDistance = windowSize;
			boolean firstMatchFound = false;
			// search for the position, not iterate everyone
			for(int j=firstMatchIndex;j<allMotifs.size();j++){
				Point motif =allMotifs.get(j);
				if (peak.getChrom().equalsIgnoreCase(motif.getChrom())){
					int distance = peak.distance(motif);
					if (distance<=nearestDistance){
						if (!firstMatchFound){		// update on the first match
							firstMatchIndex = j;
							firstMatchFound = true;
						}							
						nearestDistance = distance;
						nearestIndex = j;
					}
					else{		// if distance > nearest, and found a match already, it is getting farther away, stop search 
						if (firstMatchFound){	
							break;
						}
					}					
				}
			}
			
			if (nearestIndex !=-1){			// motif hit is within the window
				Point nearestMotif = allMotifs.get(nearestIndex);
				peaksOffsets.put(peak, peak.offset(nearestMotif));
			}
		}
		return peaksOffsets;
	}	
		
	private class TrueHit{
		Point trueHit;
		Point peak;
		int rank;	// rank of the peak in the prediction list, method dependent
		int offset;	// positive: peak is on 5' direction of motif, negative: otherwise
		TrueHit(Point trueHit, Point peak, int rank){
			this.trueHit = trueHit;
			this.peak = peak;
			this.rank = rank;
			offset = peak.offset(trueHit);
		}
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
					double leftScore= profiler.getHigherScore(windowSize-z);
					double rightScore= profiler.getHigherScore(windowSize+z);				
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
		CommonUtils.writeFile("CTCF_motif_distance_byScore_"+windowSize+"_"+method+".txt", sb.toString());
	}	
		
	private void motifBasedAnalysis(){
		long tic = System.currentTimeMillis();
		// Get the set of motif matches from all methods
		System.out.print("\nGetting motifs from GPS result ...\n");
		maps_gps = getAllNearestMotifs(peaks_gps, motifThresholds[0]);
		System.out.print("\nGetting motifs from MACS result ...\n");
		maps_macs = getAllNearestMotifs(peaks_macs, motifThresholds[0]);
		System.out.println(CommonUtils.timeElapsed(tic));
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
				motifs_gps.add(motif);
			for (Point motif: maps_macs.keySet())
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
					for (TrueHit hit:maps_gps.get(motif))
						if (hit.rank<=(j+1)*step){			// if the peak is in top rank range
							total+= Math.abs(hit.offset);
							count++;
						}
					for (TrueHit hit:maps_macs.get(motif))
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
		CommonUtils.writeFile("CTCF_GPS_SpatialResolution_"+windowSize+".txt", sb.toString());

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
		CommonUtils.writeFile("CTCF_MACS_SpatialResolution_"+windowSize+".txt", sb.toString());
		System.out.println(CommonUtils.timeElapsed(tic));
	}

	// get the all nearest motif hits in the region (within windowSize) 
	// around each peak in the list
	// if a motif have more than one peaks nearby, associate it with all the peaks
	
	// returns a map of motif -- map to all associated peaks
	private HashMap<Point, ArrayList<TrueHit>> getAllNearestMotifs(ArrayList<Point> peaks, double threshold){	
		HashMap<Point, ArrayList<TrueHit>> motifs = new HashMap<Point, ArrayList<TrueHit>>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(motif);
		
		for(int i=0; i<peaks.size(); i++){
			if (i % 1000==0)
				System.out.println(i);
			Point peak = peaks.get(i);
			Region r= peak.expand(windowSize);
			WeightMatrixScoreProfile profiler = scorer.execute(r);
			//search for whole region
			for(int z=0; z<r.getWidth(); z++){		
				double score = profiler.getHigherScore(z);
				if(score >= threshold){
					Point motifPos = new Point(genome, peak.getChrom(), r.getStart()+z+profiler.getMatrix().length()/2);
					TrueHit hit = new TrueHit(motifPos, peak, i);
					//if other peaks also associate with this motif position, add to the list
					if (motifs.containsKey(motifPos))	
						motifs.get(motifPos).add(hit) ;
					else{
						ArrayList<TrueHit> hits = new ArrayList<TrueHit>();
						hits.add(hit);
						motifs.put(motifPos, hits);
					}
				}
			}
		}
		return motifs;
	}

}
