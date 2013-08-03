package edu.mit.csail.cgs.deepseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * A BindingModel defines a (probabilistic) model of read occurrences around a binding event.
 * The probability is directional (i.e. stranded). 
 * Given a signed distance from the event, the BindingModel should return 
 * a relative probability of seeing a read at that distance.
 * 
 * @author shaunmahony
 *
 */
public class BindingModel {
	public final static int SMOOTHING_STEPSIZE = 30;
	public final static int SMOOTHING_AVG_PTS = 30;
	protected int min, max;		// the start and end position
	protected int summit;		// the position of highest prob point
	protected double[] data;
	protected double[] probs;
	protected double[] logProbs;
	private String fileName;

	protected List<Pair<Integer, Double>> empiricalDistribution;
	
	public BindingModel(File f, int minDist, int maxDist){
		min=0; max=0;
		fileName = f.getName();
		try {
			empiricalDistribution = new LinkedList<Pair<Integer,Double>>(); 
			BufferedReader reader = new BufferedReader(new FileReader(f));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>=2){
	              Integer dist = new Integer(words[0]);
	              //make sure the current data point is within the specified range
	              if ((dist.intValue() >= minDist) && (dist.intValue() <= maxDist)) {
	                Pair<Integer,Double> p = new Pair<Integer,Double>(dist, new Double(words[1]));
	                if (p.cdr().doubleValue()>=0)	// should be non-negative value
	                	empiricalDistribution.add(p);
	                else {
	                	System.err.println("\nRead distribution file contains negative probability(count) value!"); 
	                	System.exit(1);
	                }
	              }
	            }
	        }
	        loadData(empiricalDistribution);
			makeProbabilities();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public BindingModel(File f) {
		this(f, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	
	public BindingModel(List<Pair<Integer, Double>> bindingDist){
		min=0; max=0;
		empiricalDistribution=bindingDist;
		loadData(bindingDist);
		makeProbabilities();
	}
	
	//Accessors
	/** return the leftmost coordinate of the read model */
	public int getMin(){return min;}
	/** return the rightmost coordinate of the read model */
	public int getMax(){return max;}
	/** return the coordinate of summit of the read model */
	public int getSummit(){return summit;}
	/** return the larger value of positive side or negative side of the read model */
	public int getRange(){return Math.max(Math.abs(min), Math.abs(max));}
	/** return the total width of read model (positive side + negative side) */
	public int getWidth(){return data.length;}
	public double[] getProbabilities(){	return  probs.clone();}
	public double[] getLogProbabilities() { return logProbs.clone();}
	public String getFileName() {
		return fileName;
	}
	public void setFileName(String fileName) {
		this.fileName = fileName;
	}	
	//Return a pair of distances corresponding to the central probability interval provided
	//Can be used to provide hit extension lengths
	public Pair<Double,Double> probIntervalDistances(double prob){
		double ends=(1-prob)/2;
		double probSum=0;
		boolean firstFound=false, secondFound=false;
		int first=min, second=max;
		for(int i=min; i<=max; i++){
			probSum+=probability(i);
			if(!firstFound && probSum>ends){
				firstFound=true;
				first=i;
			}else if(!secondFound && probSum>(1-ends)){
				secondFound=true;
				second=i;
			}
		}
		Pair<Double,Double> intervalDists = new Pair<Double, Double>(new Double(first), new Double(second));
		return intervalDists;
	}
	//Return a pair of distances corresponding to the interval in the probability landscape that is above that expected from a uniform distribution
	//Can be used to provide hit extension lengths
	public Pair<Double,Double> probIntervalAboveUniform(){
		double prob = 1/(double)(max-min);
		boolean firstFound=false, secondFound=false;
		int first=min, second=max;
		for(int i=min; i<=max; i++){
			if(!firstFound && probability(i)>prob){
				firstFound=true;
				first=i;
			}else if(i>0 && firstFound && !secondFound && probability(i)<prob){
				secondFound=true;
				second=i;
			}
		}
		Pair<Double,Double> intervalDists = new Pair<Double, Double>(new Double(first), new Double(second));
		return intervalDists;
	}
	
	//Load data
	private void loadData(List<Pair<Integer, Double>> bindingDist){
		//Assumes the list is sorted//
		
		//Find max, min values first
		for(Pair<Integer, Double> p : bindingDist){
			if(p.car()<min)
				min=p.car();
			if(p.car()>max)
				max=p.car();
		}
		//Initialize arrays
		data = new double[(max-min)+1];
		probs = new double[(max-min)+1];
		logProbs = new double[(max-min)+1];
		for(int i=0; i<=(max-min); i++){
			data[i]=0; probs[i]=0; logProbs[i] = Double.NEGATIVE_INFINITY;
		}
		
		//Populate the data array (assumes sorted)
		int last=min-1;
		for(Pair<Integer, Double> p : bindingDist){
			int index = p.car();
			double val = p.cdr();
			//if list is not properly sorted (need to make this into an exception)
			if(index-last<0){
				System.err.println("Incorrectly sorted binding read distribution data!"); 
				System.exit(1);
			}
			//if unevenly spaced, smooth linearly between values
			if(index-last>1){
				double lastVal=dataVal(last);
				double step = (val-lastVal)/(double)(index-last);
				for(int i=1; i<(index-last); i++){
					data[(last+i)-min]=lastVal+(step*(double)i);
				}
			}
			data[index-min]=val;
			
			last = p.car();
		}
	}
	
	//Set a probability landscape according to the data. 
	public void makeProbabilities(){
		double totalVal=0;
		for(int i=min; i<=max; i++){
			totalVal+=dataVal(i);
		}
		for(int i=min; i<=max; i++){
			probs[i-min] = dataVal(i)/totalVal; 
			logProbs[i-min] = Math.log(probs[i-min]);
		}
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
		
		// update empiricalDistribution with normalized probability
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>> ();
		for(int i=min; i<=max; i++){
			newDist.add(new Pair<Integer, Double>(i, probability(i)));
		}
		empiricalDistribution=newDist;
//		for(int i=min; i<=max; i++){
//			System.out.println(i+"\t"+ probs[i-min]);
//		}
	}
	
	public void smooth(int splineStepSize, int avgStepSize){
		probs=StatUtil.cubicSpline(probs, splineStepSize, avgStepSize);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
	}
	public void smoothGaussian (int kernelWidth){
		probs=StatUtil.gaussianSmoother(probs, kernelWidth);
		Pair<Double, TreeSet<Integer>> sorted = StatUtil.findMax(probs);
		summit = sorted.cdr().first()+min;
	}	
	//Look up the probability corresponding to a distance
	// Distance should be defined as (Read position - Peak position)
	public double probability(int distance){
		if(distance<min || distance>max){
			return(0.0);
		}else{
			return(probs[distance-min]);
		}
	}
	public double probability_extended(int distance){
		if(distance<min)
			return(probs[0]);
		else if (distance>max)
			return (probs[probs.length-1]);
		else
			return(probs[distance-min]);
	}	
	
	public double logProbability(int distance) {
    if(distance<min || distance>max){
      return(Double.NEGATIVE_INFINITY);
    }else{
      return(logProbs[distance-min]);
    }	  
	}
	
	//Look up the data corresponding to a distance
	public double dataVal(int distance){
		if(distance<min || distance>max){
			return(0.0);
		}else{
			return(data[distance-min]);
		}
	}
	//Return the distance from the maximum value to the zero point
	public int maxShift(){
		int shift = 0; 
		double maxVal=0;
		for(int i=0; i<max; i++){
			if(dataVal(i)>maxVal){
				shift=i; maxVal=dataVal(i);
			}
		}return(shift);
	}
	
	// expand the model by setting prob of new positions uniformly
	// with the prob of min or max (usually very small, but not 0)
	public BindingModel getExpandedModel(int minExt, int maxExt){
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>>();
		for (int i=-minExt; i<=-1; i++)
			newDist.add(new Pair<Integer, Double>(min+i, probability(min)));
		for (int i=min; i<=max; i++)
			newDist.add(new Pair<Integer, Double>(i, probability(i)));
		for (int i=1; i<=maxExt; i++)
			newDist.add(new Pair<Integer, Double>(max+i, probability(max)));
		return new BindingModel(newDist);
	}
	
	// expand the model by setting prob of new positions "linearly"
	// with the prob of min or max (usually very small, but not 0)
	public BindingModel getLinearExpandedModel(int minExt, int maxExt){
		double smallestProb = 1e-8;
		double minProb = probability(min);
		double maxProb = probability(max);
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>>();
		for (int i=-minExt; i<=-1; i++)
			newDist.add(new Pair<Integer, Double>(min+i, smallestProb+(minProb-smallestProb)*(i+minExt)/(-1+minExt)));
		for (int i=min; i<=max; i++)
			newDist.add(new Pair<Integer, Double>(i, probability(i)));
		for (int i=1; i<=maxExt; i++)
			newDist.add(new Pair<Integer, Double>(max+i, smallestProb+(maxProb-smallestProb)*(maxExt-i)/(maxExt-1)));
		return new BindingModel(newDist);
	}	
	// expand the model by setting prob of new positions "exponentially"
	// with the prob of min or max (usually very small, but not 0)
	public BindingModel getExponentialExpandedModel(int minExt, int maxExt){
		double smallestProb = 1e-8;
		double minProb = probability(min);
		double maxProb = probability(max);
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>>();
		for (int i=-minExt; i<=-1; i++)
			newDist.add(new Pair<Integer, Double>(min+i, minProb/(-i)));
		for (int i=min; i<=max; i++)
			newDist.add(new Pair<Integer, Double>(i, probability(i)));
		for (int i=1; i<=maxExt; i++)
			newDist.add(new Pair<Integer, Double>(max+i, maxProb/i));
		return new BindingModel(newDist);
	}	
	
	//Print probs to a file
	public void printToFile(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			for(int i=min; i<=max; i++){
				fout.write(i+"\t"+probability(i)+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Command-line interface to load a BindingModel from a file
	 * Example running: <br>
	 * <tt>
	 * --in oct4.shear.ext.txt --out out_oct4.txt --out_smooth out_oct4_smooth.txt
	 * </tt>
	 */
	public static void main(String[] args){
		if(Args.parseArgs(args).contains("in")){
			String infile = Args.parseString(args, "in", null);
			String outfile = Args.parseString(args, "out", "out.model");
			String outfile_smooth = Args.parseString(args, "out_smooth", "out_smooth.model");
			File pFile = new File(infile);
			if(!pFile.isFile()){
				System.err.println("Invalid file name");
				System.exit(1);
			}
	        //File loaded, make a BindingModel
	        BindingModel model = new BindingModel(pFile);
	        model.printToFile(outfile);
	        
	        model.smooth(SMOOTHING_STEPSIZE, SMOOTHING_AVG_PTS);
	        model.printToFile(outfile_smooth);
		}else{
			System.out.println("Usage: BindingModel --in GPSfileName --out outfile");
		}
	}
	
	public List<Pair<Integer, Double>> getEmpiricalDistribution() {
		List<Pair<Integer, Double>> newDist = new ArrayList<Pair<Integer, Double>> ();
		for (Pair<Integer, Double> p: empiricalDistribution)
			newDist.add(p);
		return newDist;
	}
	
	public int findNewMax(){
		double leftProb = probability(min);
		double rightProb = probability(max);
		// if left side is higher
		if ((leftProb-rightProb)/rightProb > 0.5){
			for (int i=probs.length-1;i>=0;i--){
				if (probs[i]>leftProb){
					return max-(probs.length-1-i)/2;
				}
			}
		}
		// if right side is higher
		if ((rightProb-leftProb)/leftProb > 0.5){
			for (int i=0; i<probs.length;i++){
				if (probs[i]>rightProb){
					return max+i/2;
				}
			}
		}
		return max;
	}

	/**
	 * Estimate a better read profile ranges
	 * @param minLeft minimum range on the left side
	 * @param minRight	minimum range on the right side
	 * @return The estimated profile ranges <left, right>
	 */
	public Pair<Integer, Integer> getNewEnds(int minLeft, int minRight){
		double halfHeight = probability(summit)*0.5;
		int leftHalfHeightEnd = 0;
		int rightHalfHeightEnd = 0;
		for (int i=min; i<=max; i++){
			if (probability(i)>=halfHeight){
				leftHalfHeightEnd = i;
				break;
			}
		}
		for (int i=max; i<=min; i--){
			if (probability(i)>=halfHeight){
				rightHalfHeightEnd = i;
				break;
			}
		}
		int left=Math.max(minLeft, Math.abs(summit-(summit-leftHalfHeightEnd)*4));
		int right=Math.max(minRight, Math.abs(summit+(rightHalfHeightEnd-summit)*3));

		return new Pair<Integer, Integer>(left, right);
	}
	// shift two array elements to give best KL-divergence
	// it will mutate the two input arrays.
	// assume a and b have same length
	public static int minKL_Shift (double[] a, double[] b){
		return minKL_Shift (a,b,20);
	}
	/** Shift two array elements to give best KL-divergence
	 * it will mutate the two input arrays.
	 * assume a and b have same length
	 * 
	 * @param a
	 * @param b
	 * @param range the range to shift
	 * @return
	 */
	public static int minKL_Shift (double[] a, double[] b, int range){
		if (a.length!=b.length)
			return 0;
		int length = a.length;
		double minKL=Double.MAX_VALUE;
		int shift=0;
		double[] a_new = new double[length];
		double[] b_new = new double[length];
		double[] a1 = new double[length];
		double[] b1 = new double[length];
		
		for (int i=-range;i<range;i++){
			if (i<0){				//a shift forward, b shift backward
				System.arraycopy(a, 0, a1, -i, length+i);
				System.arraycopy(a, length+i, a1, 0, -i);
				System.arraycopy(b, -i, b1, 0, length+i);
				System.arraycopy(b, 0, b1, length+i, -i);
			}
			else if (i>0){ 			//a shift backward, b shift forward
				System.arraycopy(a, i, a1, 0, length-i);
				System.arraycopy(a, 0, a1, length-i, i);
				System.arraycopy(b, 0, b1, i, length-i);
				System.arraycopy(b, length-i, b1, 0, i);
			}
			else{
				System.arraycopy(a, 0, a1, 0, length);
				System.arraycopy(b, 0, b1, 0, length);
			}

			double kl = StatUtil.log_KL_Divergence(a1, b1);

			if (kl<minKL){
				minKL=kl;
				shift=i;
				System.arraycopy(a1, 0, a_new, 0, length);
				System.arraycopy(b1, 0, b_new, 0, length);
//				System.out.println("KL= "+minKL+" when shifting "+shift+" bp");
			}
		}

		System.arraycopy(a_new, 0, a, 0, length);
		System.arraycopy(b_new, 0, b, 0, length);

		return shift;
	}
}
