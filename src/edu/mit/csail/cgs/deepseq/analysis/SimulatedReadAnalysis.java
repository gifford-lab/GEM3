package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.StrandedBase;
import edu.mit.csail.cgs.deepseq.discovery.BindingMixture;
import edu.mit.csail.cgs.deepseq.features.ComponentFeature;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.deepseq.utilities.ReadSimulator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.stats.VectorUtil;

public class SimulatedReadAnalysis {
	
	private String[] args;
	private int numSamples;
	private int randSeed;
	private int noiseRandSeed;
	private double noise;
	private int analysisType;
	
	private ReadSimulator sim;
	private BindingModel model;
	private BindingMixture mixModel;
	private String timeString;
	private double[] gaussian;
	private int smoother_width=25;
	
	private Vector<Double> resolutions = new Vector<Double>();
	private Vector<Double> logKLs = new Vector<Double>();
	private double[] meanResolutions;
	private double[] stdevResolutions;
	private double[] meanLogKLs;
	private double[] stdevLogKLs;

	private int[] readCounts;
	private int[] readCounts0 = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
	private int[] readCounts1 = {10, 20, 30, 40, 70, 100, 150, 200, 300, 400, 500, 700, 1000};
	private int[] readCounts2 = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200};
	private int[] readCounts3 = {100};
	/**
	 * @param args
	 * 
	 * Example of parameters
--read_distribution "Read_Distribution_default.txt"
--geninfo "sim.geninfo.txt" 
--analysisType 0
--readCountOption 1
--numSamples 100
--randSeed 0
--noiseRandSeed 1
--noise 0
	 */
	public static void main(String[] args) {
		SimulatedReadAnalysis analysis = new SimulatedReadAnalysis(args);
		int type = Args.parseInteger(args, "analysisType", 0);
		switch (type){
		case 0:	// mixture model prediction, spatial resolution and peak shape
			analysis.GPS_RobustnessAnalysis();
			analysis.output("GPS_Robustness");
			break;
		case 1:	// get logKL distribution with different value of Gaussian Kernel width
//			int[] width = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 30};
			int[] width = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
			for (int i:width){
				analysis.setGaussianWidth(i);
				analysis.run_shape_Analysis();
				analysis.output("Shape_width");
			}
			break;
		case 2:	// print all the logKL data for all samplings
			analysis.peakShapeData();
			break;
		case 3: // get peak shape logKL distribution parameter, mean and std (KL ~ log-nomal distribution)
			analysis.peakShapeDistributionParameters();
			analysis.output("PeakShapeDistribution");
			break;
    	default: 
    		System.err.println("Unknown analysis type "+type);
		}
		

	}
	 
	public SimulatedReadAnalysis(String[] args){ 
		this.args = args;
		
		String modelFile = Args.parseString(args, "read_distribution", null);

		
		for (String arg: args){
			System.out.print(arg+" ");
		}
		System.out.println();
		
		File mFile = new File(modelFile);
		if(!mFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
		
		analysisType = Args.parseInteger(args, "analysisType", 0);
		numSamples = Args.parseInteger(args, "numSamples", 100);
		randSeed = Args.parseInteger(args, "randSeed", 0);
		noiseRandSeed = Args.parseInteger(args, "noiseRandSeed", 1);
		noise = Args.parseDouble(args, "noise", 0.0);
		
		int option = Args.parseInteger(args, "readCountOption", 0);
		switch (option){
		case 0:	
			readCounts = readCounts0;
			break;
		case 1:
			readCounts = readCounts1;
			break;
		case 2:	
			readCounts = readCounts2;
			break;
		case 3: 
			readCounts = readCounts3;
			break;
    	default: 
    		System.err.println("Unknown readCountOption "+option);
		}
        //File loaded, make a BindingModel
        model = new BindingModel(mFile);
        
        //init the Guassian kernel prob. for smoothing the read profile of called events
		gaussian = new double[model.getWidth()]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, smoother_width*smoother_width);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);

        //Initialize the ReadSimulator
        sim = new ReadSimulator(model, Args.parseString(args, "geninfo", ""));
        sim.setRandSeed(randSeed);		// set seed number for random generator
        sim.setNoiseRandSeed(noiseRandSeed);
		sim.setNoiseProb(noise);        

		meanResolutions = new double[readCounts.length];
		stdevResolutions = new double[readCounts.length];
		meanLogKLs = new double[readCounts.length];
		stdevLogKLs = new double[readCounts.length];
	}
	
	// Simulated reads
	// predict peaks using GPS
	// Get stats: resulutions, logKLs
	public void GPS_RobustnessAnalysis()
	{
        mixModel = new BindingMixture(args, true);
		long tic = System.currentTimeMillis();

		for (int j=0;j<readCounts.length;j++){
			if (readCounts[j]==0)
				continue;
	        for (int i=0;i<numSamples; i++){
	        	sim.restart();
	        	List<ReadHit> reads = sim.simulateBothStrands(readCounts[j]);
	        	
	        	ArrayList<ComponentFeature> peaks = 
	        		mixModel.simpleExecute(StrandedBase.convertReadHits(reads), sim.getSimRegion());
//	        	System.out.print(peaks.size());
	        	if (peaks.size()!=1)
	            	System.err.println("Found " +peaks.size()+ " events.");
//	        		{System.err.print("\n"+peaks.size());
//	        	}
	        	for (ComponentFeature peak: peaks){
	        		double result = peak.getPeak().getLocation();
	        		resolutions.add(Math.abs(result-sim.DEFAULT_PEAK_LOCATION));
	        		logKLs.add(peak.getAverageLogKL());
//	        		System.out.println(result+"\t"+(result-sim.DEFAULT_START));
	        	}
	        }
	        if (!resolutions.isEmpty()){
		        meanResolutions[j]=VectorUtil.mean(resolutions);
		        stdevResolutions[j]=VectorUtil.stdev(resolutions);
	        }
	        if (!logKLs.isEmpty()){
		        meanLogKLs[j]=VectorUtil.mean(logKLs);
		        stdevLogKLs[j]=VectorUtil.stdev(logKLs);
	        }
	        resolutions.removeAllElements();
	        logKLs.removeAllElements();
		}
		timeString = "Total time:\t"+CommonUtils.timeElapsed(tic);
        System.out.println(timeString);
	}

	public void run_shape_Analysis()
	{
        mixModel = new BindingMixture( args, false);
		long tic = System.currentTimeMillis();

		for (int j=0;j<readCounts.length;j++){
			if (readCounts[j]==0)
				continue;
	        for (int i=0;i<numSamples; i++){
	        	sim.restart();
	        	List<ReadHit> reads = sim.simulateBothStrands(readCounts[j]);
	        	
	        	// select what type of analysis to carry on
		        switch(analysisType){
		        	case 1: 
		        	case 3:
		        		peakShapeAnalysis(reads);break;
		        	default: 
		        		System.err.println("Unknown analysis type "+analysisType);
	        	}
	        }
	        if (!resolutions.isEmpty()){
		        meanResolutions[j]=VectorUtil.mean(resolutions);
		        stdevResolutions[j]=VectorUtil.stdev(resolutions);
	        }
	        if (!logKLs.isEmpty()){
		        meanLogKLs[j]=VectorUtil.mean(logKLs);
		        stdevLogKLs[j]=VectorUtil.stdev(logKLs);
	        }
	        resolutions.removeAllElements();
	        logKLs.removeAllElements();
		}
		timeString = "Total time:\t"+CommonUtils.timeElapsed(tic);
        System.out.println(timeString);
	}
	
	public void output(String name_prefix){
		CommonUtils.printArray(readCounts, "", "\n");
		CommonUtils.printArray(meanResolutions, "", "\n");
		CommonUtils.printArray(stdevResolutions, "", "\n");
		CommonUtils.printArray(meanLogKLs, "", "\n");
		CommonUtils.printArray(stdevLogKLs, "", "\n");
        
        StringBuilder output = new StringBuilder();
        output.append("#args: ");
        for (String arg: args)
			output.append(arg+" ");
		output.append("\n#");
//		output.append(mixModel.getConfigString()).append("\n");
		output.append(timeString).append("\n");
		output.append("Gaussian width:\t").append(smoother_width).append("\n");
		output.append("Sample number:\t").append(numSamples).append("\n\n");
		
        output.append(CommonUtils.arrayToString(readCounts)).append("\n")
        	  .append(CommonUtils.arrayToString(meanResolutions)).append("\n")
        	  .append(CommonUtils.arrayToString(stdevResolutions)).append("\n")
        	  .append(CommonUtils.arrayToString(meanLogKLs)).append("\n")
        	  .append(CommonUtils.arrayToString(stdevLogKLs));
        
       CommonUtils.writeFile(name_prefix + "_noise"+String.format("%.3f", noise)+ "_width"+smoother_width+".txt", output.toString());
	}
	

	private void setGaussianWidth(int width){
		smoother_width = width;
		NormalDistribution gaussianDist = new NormalDistribution(0, width*width);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
	}
	
	// estimate KL peak shape parameter using a group of reads
	private void peakShapeAnalysis(List<ReadHit> reads){
		double logKL_plus, logKL_minus;
		// array plus and minus to store read counts along the range of binding model
		double[] plus = new double[model.getWidth()];
		double[] minus = new double[model.getWidth()];
		for (ReadHit read:reads){
			if (read.getStrand()=='+'){
				int index = read.getFivePrime()-sim.DEFAULT_PEAK_LOCATION-model.getMin();
				if (index<0 || index>minus.length-1)
					System.out.print("");
				plus[index]++;
				continue;
			}
			if (read.getStrand()=='-'){
				int index = sim.DEFAULT_PEAK_LOCATION-read.getFivePrime()-model.getMin();
				if (index<0 || index>minus.length-1)
					System.out.print("");
				minus[index]++;
			}
		}
		logKL_plus = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(plus, gaussian));
		logKL_minus = StatUtil.log_KL_Divergence(model.getProbabilities(), StatUtil.symmetricKernelSmoother(minus, gaussian));
//		resolutions.add(Math.abs(logKL_plus - logKL_minus)/(logKL_plus + logKL_minus));	// use resolutions list to store symmetric score
		logKLs.add((logKL_plus+logKL_minus)/2);
	}
	
	private void peakShapeData(){
		int readNum = 20;
        for (int i=0;i<numSamples; i++){
        	sim.restart();
        	List<ReadHit> reads = sim.simulateBothStrands(readNum);        	
        	peakShapeAnalysis(reads);
        }
        StringBuilder output = new StringBuilder();
        for (double kl: logKLs){
        	output.append(kl).append("\n");
        }
        CommonUtils.writeFile("Read"+readNum+"_logKLs.txt", output.toString());
	}
	
	private void peakShapeDistributionParameters(){
		int count = 100;
		readCounts = new int[count];
		meanResolutions = new double[readCounts.length];
		stdevResolutions = new double[readCounts.length];
		meanLogKLs = new double[readCounts.length];
		stdevLogKLs = new double[readCounts.length];
		
		for (int r=5;r<=count;r++){
			readCounts[r-1]=r;
		}
		
		this.run_shape_Analysis();
	}
}
