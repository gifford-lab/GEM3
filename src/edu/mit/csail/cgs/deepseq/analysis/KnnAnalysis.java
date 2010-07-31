package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSParser;
import edu.mit.csail.cgs.ewok.verbs.chipseq.GPSPeak;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class KnnAnalysis {
	
	private int kk;	
	private double fdr_threshold;

	private Genome genome;
	private String[] args;
	protected String outName="out";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		KnnAnalysis analysis = new KnnAnalysis(args);
		analysis.run();
	}
	
	KnnAnalysis(String[] args){
		this.args = args;

		try {
			genome = Organism.findGenome("mm8");
		} catch (NotFoundException e) {
			e.printStackTrace();
		}

		// some parameters
		kk = Args.parseInteger(args, "K", 100)+1;	// add 1 for counting self as nearest
		fdr_threshold = Args.parseDouble(args, "fdr", 0.015);
		outName = Args.parseString(args,"out",outName);
	}

	
	private void run(){
		List<GPSPeak> ipPeaks=null;
		ArrayList<KnnPoint> knnPoints=new ArrayList<KnnPoint>();
		List<GPSPeak> ctrlPeaks=null;
		
		// load GPS results from IP and Control
		String filename = Args.parseString(args, "IP", null);
		if (filename==null){
			System.err.println("GPS IP file is not found.");
			System.exit(-1);
		}
		File gpsFile = new File(filename);
		ipPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
		for (GPSPeak p: ipPeaks){
			knnPoints.add(new KnnPoint(p, true));
		}	
		
		filename = Args.parseString(args, "Ctrl", null);
		if (filename==null){
			System.err.println("GPS Control file is not found.");
			System.exit(-1);
		}
		gpsFile = new File(filename);
		ctrlPeaks = GPSParser.parseGPSOutput(gpsFile.getAbsolutePath(), genome);
		for (GPSPeak p: ctrlPeaks){
			knnPoints.add(new KnnPoint(p, false));
		}
		
		// variance-normalize data
		double strength[]=new double[knnPoints.size()];
		double shape[]=new double[knnPoints.size()];
		double fold[]=new double[knnPoints.size()];
		for (int i=0;i<knnPoints.size();i++){
			KnnPoint p = knnPoints.get(i);
			strength[i]=p.strength;
			shape[i]=p.shape;
			fold[i]=p.fold;
		}
		double varStrength = StatUtil.std(strength);
		double varShape = StatUtil.std(shape);
		double varFold = StatUtil.std(fold);
		for (int i=0;i<knnPoints.size();i++){
			KnnPoint p = knnPoints.get(i);
			p.strength=p.strength/varStrength;
			p.shape=p.shape/varShape;
			p.fold=p.fold/varFold;
		}		
		
		// distance calculation
		long tic = System.currentTimeMillis();
		calcDistance_Naive(knnPoints);
		System.out.println(CommonUtils.timeElapsed(tic));
		

//		int[] kValues = new int[10];
//		for (int i=0; i<10; i++){
//			kValues[i] = kk*(10-i)/10;
//		}

		int[] kValues = new int[1];
		kValues[0] = kk-1;
		
		// output counts for FDR calculation
		double fdr_cutoff=0;
		int n_cutoff=0;
		for (int j=0;j<kValues.length;j++){
			int k = kValues[j];
			double[] truePositives = new double[k+1];
			double[] falsePositives = new double[k+1];
			for (int n=1;n<=k;n++){
				double truePositive=0;	// ip points that are called
				double falsePositive=0; // ctrl points that are called
				for (KnnPoint p: knnPoints){
					if (p.isIP)
						if (p.getIPCount(k)>=n)
							truePositive++;
					if (!p.isIP)
						if (p.getIPCount(k)>=n)
							falsePositive++;
				}
				truePositives[n]=truePositive;
				falsePositives[n]=falsePositive;
				double fdr = falsePositive/truePositive;
				if (fdr_threshold >= fdr){
					if (fdr_cutoff < fdr){
						fdr_cutoff = fdr;
						n_cutoff = n;
					}
				}
			}
			System.out.println();
			StringBuilder sb = new StringBuilder();
			for (int n=1;n<=k;n++){
				sb.append(n).append("\t").append(truePositives[n]).append("\t")
					.append(falsePositives[n]).append("\t")
					.append(String.format("%.3f", falsePositives[n]/truePositives[n]))
					.append("\n");
			}
			CommonUtils.writeFile(outName+"_K"+k+".txt", sb.toString());
			
			// output n-sorted peaks list
			Collections.sort(knnPoints, new Comparator<KnnPoint>(){
			    public int compare(KnnPoint o1, KnnPoint o2) {
			        return o1.compareByNandP(o2);
			    }
			});
			sb = new StringBuilder();
			for (KnnPoint p:knnPoints){
				if(p.isIP && p.n>=n_cutoff)
					sb.append(p.peak.toGPS()).append("\t")
					.append(p.n).append("\t").append(k).append("\n");
			}
			CommonUtils.writeFile(outName+"_K"+k+"_N"+n_cutoff+"_Peaks.txt", sb.toString());
			sb = new StringBuilder();
			for (KnnPoint p:knnPoints){
				if(p.isIP && p.n<n_cutoff)
					sb.append(p.peak.toGPS()).append("\t")
					.append(p.n).append("\t").append(k).append("\n");
			}
			CommonUtils.writeFile(outName+"_K"+k+"_N"+n_cutoff+"_Insig_Peaks.txt", sb.toString());
		}
	}
	// naive implementation: calculating distance to all other points
	private void calcDistance_Naive(ArrayList<KnnPoint> knnPoints){
		System.out.println("Total size: "+knnPoints.size());
		ArrayList<KnnPoint> tmpPoints=new ArrayList<KnnPoint>();
		int count = 1;
		tmpPoints.addAll(knnPoints);
		for (KnnPoint p: knnPoints){
			if (count % 1000 == 0) {
			    if (count % 10000 == 0) {
			      System.out.println(count);
			    }
			    else {
			      System.out.print(count);
			    }
			}
			else if (count % 500 == 0) {
			    System.out.print(".");
			}
			for (KnnPoint q: tmpPoints){
				q.calcDistance(p);
			}
			Collections.sort(tmpPoints);
			boolean[] isIP = new boolean[kk];
			for (int i=0;i<kk;i++){
				isIP[i]=tmpPoints.get(i).isIP;
			}
			p.setNeighbors(isIP);
			count++;
		}
	}
	
	// the basic idea is to avoid calculating distance to all other points, 
	// but only calculate distance for the intersecting bands around the point
	private void calcDistance_BandIntersect(ArrayList<KnnPoint> knnPoints){
		ArrayList<KnnPoint> strengthSorted=new ArrayList<KnnPoint>();
		strengthSorted.addAll(knnPoints);
		Collections.sort(strengthSorted, new Comparator<KnnPoint>(){
		    public int compare(KnnPoint o1, KnnPoint o2) {
		        return o1.compareByStrength(o2);
		    }
		});
		double[] strengthVector = new double[strengthSorted.size()];
		for (int i=0; i<strengthSorted.size(); i++){
			strengthVector[i]=strengthSorted.get(i).strength;
		}
		ArrayList<KnnPoint> shapeSorted=new ArrayList<KnnPoint>();
		shapeSorted.addAll(knnPoints);
		Collections.sort(shapeSorted, new Comparator<KnnPoint>(){
		    public int compare(KnnPoint o1, KnnPoint o2) {
		        return o1.compareByShape(o2);
		    }
		});
		double[] shapeVector = new double[shapeSorted.size()];
		for (int i=0; i<shapeSorted.size(); i++){
			shapeVector[i]=shapeSorted.get(i).shape;
		}
		ArrayList<KnnPoint> foldSorted=new ArrayList<KnnPoint>();
		foldSorted.addAll(knnPoints);
		Collections.sort(foldSorted, new Comparator<KnnPoint>(){
		    public int compare(KnnPoint o1, KnnPoint o2) {
		        return o1.compareByFold(o2);
		    }
		});
		double[] foldVector = new double[foldSorted.size()];
		for (int i=0; i<foldSorted.size(); i++){
			foldVector[i]=foldSorted.get(i).fold;
		}

		for (KnnPoint p:knnPoints){
			int overlap = 0;
			int bandwidth = kk*4;
			Set<KnnPoint> union=null;
			while(overlap<kk){
				int[] edges_strength = findEdges(strengthVector, p.strength, bandwidth);
				Set<KnnPoint> points_strength= new HashSet<KnnPoint>();
				for (int i=edges_strength[0]; i<=edges_strength[1]; i++)
					points_strength.add(strengthSorted.get(i));
				int[] edges_shape = findEdges(shapeVector, p.shape, bandwidth);
				Set<KnnPoint> points_shape= new HashSet<KnnPoint>();
				for (int i=edges_shape[0]; i<=edges_shape[1]; i++)
					points_shape.add(shapeSorted.get(i));
				int[] edges_fold = findEdges(foldVector, p.fold, bandwidth);
				Set<KnnPoint> points_fold= new HashSet<KnnPoint>();
				for (int i=edges_fold[0]; i<=edges_fold[1]; i++)
					points_fold.add(foldSorted.get(i));
				SetTools<KnnPoint> tool = new SetTools<KnnPoint>();
				overlap = tool.intersection(points_fold, tool.intersection(points_strength, points_shape)).size();
				if (overlap>=kk){
					points_strength.clear();
					for (int i=edges_strength[2]; i<=edges_strength[3]; i++)
						points_strength.add(strengthSorted.get(i));
					points_shape.clear();
					for (int i=edges_shape[2]; i<=edges_shape[3]; i++)
						points_shape.add(shapeSorted.get(i));
					points_fold.clear();
					for (int i=edges_fold[2]; i<=edges_fold[3]; i++)
						points_fold.add(foldSorted.get(i));
					union = tool.union(points_fold, tool.union(points_strength, points_shape));
				}
				else
					bandwidth +=2*kk;
			}
			//System.out.print(bandwidth+" ");
			ArrayList<KnnPoint> tmpPoints=new ArrayList<KnnPoint>();
			tmpPoints.addAll(union);
			for (KnnPoint q: tmpPoints){
				q.calcDistance(p);
			}
			Collections.sort(tmpPoints);
			
			//****************************
//			ArrayList<KnnPoint> tmpPoints2=new ArrayList<KnnPoint>();
//			tmpPoints2.addAll(knnPoints);
//			for (KnnPoint q: tmpPoints2){
//				q.calcDistance(p);
//			}
//			Collections.sort(tmpPoints2);
//			
//			for (int i=0;i<kk;i++){
//				double distance = tmpPoints.get(i).distance;
//				double distance2 = tmpPoints2.get(i).distance;
//				
//				if (distance != distance2){
//					System.out.println(i+"\t"+distance+"\t"+distance2);
//				}
//			}
			
			//********************************
			
			
			boolean[] isIP = new boolean[kk];
			for (int i=0;i<kk;i++){
				isIP[i]=tmpPoints.get(i).isIP;
			}
			p.setNeighbors(isIP);
		}
	}
	
	// fine the edges of a band around "middle" point, as index of the "array"
	// the width of the band should be larger or equal to "width"
	// the array is assumed to be sorted in ascending order
	// the bandwidth is not exactly tight, but the whole intersect method is an approximate, so no need to be exact here
	private int[] findEdges(double[] array, double center, int width){
		// edges: the indexes to return, first 2 are used for search overlaps
		// 		  the last 2 indexes are more conservative range for calculation
		int[] edges = new int[4];
		int lowIdx=0;
		int highIdx=0;
		int idx = Arrays.binarySearch(array, center);
		lowIdx = Math.max(0, idx-width/2);
		if (lowIdx==0){
			edges[0]=0;
			edges[1]=width;
			edges[2]=0;
			edges[3]=width;
			return edges;
		}
		highIdx = Math.min(array.length-1, idx+width/2);
		if (highIdx==array.length-1){
			edges[0]=highIdx-width;
			edges[1]=highIdx;
			edges[2]=highIdx-width;
			edges[3]=highIdx;
			return edges;
		}
		if (center-array[lowIdx]>array[highIdx]-center){
			// if low edge is wider, adjust high edge to have the same width(value)
			edges[0]=lowIdx;
			highIdx = Arrays.binarySearch(array, center-array[lowIdx]+center);
			if (highIdx<0)
				highIdx = -highIdx;
			edges[1]=highIdx;
			// find another edge that will be at least as far
			edges[3]=highIdx;
			lowIdx = Arrays.binarySearch(array, center+center-array[highIdx]);
			if (lowIdx<0)
				lowIdx = -lowIdx-2;
			else
				lowIdx-=2;
			lowIdx = Math.max(0, lowIdx);
			edges[2]=lowIdx;
		}else{
			edges[1]=highIdx;
			lowIdx = Arrays.binarySearch(array, center-array[highIdx]+center);
			if (lowIdx<0)
				lowIdx = -lowIdx-2;
			else
				lowIdx-=2;
			lowIdx = Math.max(0, lowIdx);
			edges[0]=lowIdx;
			edges[2]=lowIdx;
			highIdx = Arrays.binarySearch(array, center-array[lowIdx]+center);
			if (highIdx<0)
				highIdx = -highIdx;
			edges[3]=highIdx;
		}			
		return edges;
	}
	
	class KnnPoint implements Comparable<KnnPoint>{
		GPSPeak peak; 
		boolean isIP;
		double strength;
		double shape;
		double distance;	// distance to current point
		double fold;		// IP/Ctrl
		boolean[] neighbors;
		int n;				// number of neighbors in K points
		
		KnnPoint( GPSPeak peak, boolean isIP){
			this.peak = peak;
			this.isIP = isIP;
			this.strength = peak.getStrength();
			this.shape = peak.getShape();
			if (isIP){
				fold = (peak.getStrength()+1)/(peak.getControlStrength()+1);
			}
			else{
				fold = 1;
			}
		}
		void calcDistance(KnnPoint p){
			distance = (shape-p.shape)*(shape-p.shape)+(fold-p.fold)*(fold-p.fold)+
						(strength-p.strength)*(strength-p.strength);
		}
		void setNeighbors(boolean[]isIP){
			neighbors = isIP;
		}
		// get number of neighbors in k points
		int getIPCount(int k){
			int count=0;
			k = Math.min(k, neighbors.length-1);	// length of neighbor vector should be k+1
			int start=0; 
			int end=k-1;
			if (isIP){				// excluding this IP point itself, start from 1
				start=1; end=k;
			}
			for (int i=start;i<=end;i++){
				if (neighbors[i])
					count++;
			}
			n=count;
			return count;
		}

		//Comparable default method
		public int compareTo(KnnPoint p) {	// ascending
			double diff = distance-p.distance;
			return diff==0?0:(diff<0)?-1:1;
		}
		public int compareByShape(KnnPoint p) {// ascending
			double diff = shape-p.shape;
			return diff==0?0:(diff<0)?-1:1;
		}		
		public int compareByStrength(KnnPoint p){// ascending, for finding range
			double diff = strength-p.strength;
			return diff==0?0:(diff<0)?-1:1;
		}
		public int compareByFold(KnnPoint p){// ascending, for finding range
			double diff = fold-p.fold;
			return diff==0?0:(diff<0)?-1:1;
		}
		// first sorted by number of IP neighbors
		// if tie, then sorted by p-value of the peak
		public int compareByNandP(KnnPoint p){
			double diff = n-p.n;
			if (diff==0){	
				double diff2 = peak.getPvalue()- p.peak.getPvalue();
				return	diff2==0?0:(diff2<0)?1:-1;	//p-value: descending
			}
			return diff<0?1:-1;						//n: descending
		}
	}
}
