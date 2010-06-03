package edu.mit.csail.cgs.deepseq.discovery;

import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.BindingModelGenerator;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ChipSeqPeakFinder extends StatisticalPeakFinder{

	public ChipSeqPeakFinder(DeepSeqExpt signal) {
		super(signal, null);
	}
	public ChipSeqPeakFinder(DeepSeqExpt signal, DeepSeqExpt ctrl) {
		super(signal, ctrl);
	}
	public ChipSeqPeakFinder(String[] args){
		super(args);
	}
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		
		System.out.println("Welcome to the ChipSeqFinder!!\n\n");
		
		long start = System.currentTimeMillis();
		ChipSeqPeakFinder finder = new ChipSeqPeakFinder(args);
		finder.setStrandedFinding(false);
		System.out.println("Finding enriched regions");
		List<Feature> peaks = finder.execute(); 
			
		finder.printFeatures();
		finder.printFeatures(false);
			
		if(Args.parseArgs(args).contains("warpcoords")){
			finder.printWarpCoords();
			finder.printWarpCoords(false);
		}if(Args.parseArgs(args).contains("printgff")){
			finder.printGFF();
			finder.printGFF(false);
		}
		if(Args.parseArgs(args).contains("outByLocation"))
		{
			String outNameByLocation = Args.parseString(args, "outByLocation", "");
			if(!outNameByLocation.equals(""))
			{
				finder.setOutName(outNameByLocation);
				//finder.orderFeatsByLocation();
				finder.printFeatures();
				finder.printFeatures(false);					
			}
			else
			{
				System.err.println("You have to provide explicitly a prefix name " +
						"for the peak file which will be ordered w.r.t. peaks location");
			}
		}
			
			
		System.out.println("Output file root: "+finder.getOutName());
		if(!Args.parseArgs(args).contains("noseqs")){
			finder.printPeakSequences();
		}
		long stop = System.currentTimeMillis();
		long duration = (stop - start)/1000;
		long minutes = duration/60;
		System.out.println("It took " + minutes + " minutes, " + (duration-60*minutes) + " seconds to run the algorithm.\n\n");
		
		finder.cleanup();
	}//end of main method
	
	private static void multipleRuns() {
		String[] nameExp = {
				"D0_PPARg_Day0_1",
				"D1_PPARg_Day1_1",
				"D2_PPARg_Day2_1",
				"D3_PPARg_Day3_1",
				"D4_PPARg_Day4_1",
				"D6_PPARg_Day6_1",
				"D0_RXR_Day0_1",
				"D1_RXR_Day1_1",
				"D2_RXR_Day2_1",
				"D3_RXR_Day3_1",
				"D4_RXR_Day4_1",
				"D6_RXR_Day6_1",
		};
		
		String[][] args_arr =
		{
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[0] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[0], "--noShowGeneAnnotations", "--noShowOtherAnnotations", "--outByLocation", "outLoc", "--buildmodel"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[1] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[1], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[2] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[2], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[3] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[3], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[4] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[4], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[5] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[5], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[6] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[6], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[7] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[7], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[8] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[8], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[9] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[9], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[10] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[10], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
				{"--species", "Mus musculus;mm8", "--dynback", "5000", "--dynback", "10000", "--readlen", "32", "--model", "oct4.shear.ext.txt", "--modelpeak",  "--expt", "./data/RXR/HS_Adipogen_" + nameExp[11] + ".bowtie.align", "--format", "BOWTIE", "--out", "gio_" + nameExp[11], "--noShowGeneAnnotations", "--noShowOtherAnnotations"},
		};
			
		System.out.println("Welcome to the ChipSeqFinder!!\n\n");
		
		for(int i = 0; i < args_arr.length; i++)
		{
			String[] argss = args_arr[i];
			
			long start = System.currentTimeMillis();
			System.out.println("Loading experiment: " + nameExp[i] + "...\n--------------------------------------");
			ChipSeqPeakFinder finder = new ChipSeqPeakFinder(argss);
			finder.setStrandedFinding(false);
			System.out.println("Finding enriched regions");
			List<Feature> peaks = finder.execute(); 
			
			finder.printFeatures();
			finder.printFeatures(false);
			
			
			if(Args.parseArgs(argss).contains("outByLocation"))
			{
				String outNameByLocation = Args.parseString(argss, "outByLocation", "");
				if(!outNameByLocation.equals(""))
				{
					finder.setOutName(outNameByLocation);
					//finder.orderFeatsByLocation();
					finder.printFeatures();
					finder.printFeatures(false);					
				}
				else
				{
					System.err.println("You have to provide explicitly a prefix name " +
							"for the peak file which will be ordered w.r.t. peaks location");
				}
			}
			
			
			
			long stop = System.currentTimeMillis();
			long duration = (stop - start)/1000;
			long minutes = duration/60;
			System.out.println("Exp: " + nameExp[i] + ". It took " + minutes + " minutes, " + (duration-60*minutes) + " seconds to run the algorithm.\n\n");			
		}
				
		//System.out.println("Output file root: "+finder.getOutName());
		//finder.printPeakSequences();
		
	}//end of multipleRuns
	
	public List<Feature> execute(){
		//Initialize background models
		initializeBackgrounds();
		
		//Print some more info
		System.out.println("First Round\n-----------");
        System.out.println("Sliding Window: Size="+binWidth+" Offset="+binStep);
        System.out.println("Read Length: "+readLength+", Shift: "+readShift+", 5'Ext: "+read5PrimeExt+", 3'Ext: "+read3PrimeExt);
        System.out.println("Genome-wide Poisson Threshold (Signal no less than): "+signalBacks.getGenomicModelThreshold("high"));
        if(!noControl)
        	System.out.println("Genome-wide Poisson Threshold (Control no more than): "+ctrlBacks.getGenomicModelThreshold("high"));
        System.out.println("Per-Base Poisson Threshold (Signal): "+signalPerBaseBack.getMaxThreshold('.'));
        if(!noControl)
        	System.out.println("Per-Base Poisson Threshold (Control): "+ctrlPerBaseBack.getMaxThreshold('.'));
        System.out.println("Scaling Factor: "+(noControl ? 1 : control.getScalingFactor()));
        
		//Initialize the iterator for the test regions/genes
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(gen);
		
		//Call the peak finder for the first round (finds enriched regions for the meta-peak and sets the control scaling factor
		callEnrichedRegions(testRegions, true, true);
		System.out.println("First round: "+signalPeaks.size()+" enriched signal peaks, "+controlPeaks.size()+" enriched control peaks.");
		
		//make the meta-peak
		if(buildBindingModel){
			BindingModelGenerator bmg = new BindingModelGenerator(gen, signal, signalPeaks, controlPeaks);
			bmg.setPFilter(1e-7);
			bmg.setOverRepFilter(10);
			BindingModel bm = bmg.execute();
			bm.smooth(BindingModel.SMOOTHING_STEPSIZE);
			setBindingModel(bm); //Update shift, ext
			setLRBalPeaks(false);
			setModelPeaks(true);
			bmg.printModel(new String(getOutName()+"_model.txt"));
		}
		
		//Re-initialize background models
		//setHighLogConf(tmphigh); setLowLogConf(tmplow);
		initializeBackgrounds();
		//Print some more info again
		System.out.println("Second Round\n------------");
        System.out.println("Sliding Window: Size="+binWidth+" Offset="+binStep);
        System.out.println("Read Length: "+readLength+", Shift: "+readShift+", 5'Ext: "+read5PrimeExt+", 3'Ext: "+read3PrimeExt);
        System.out.println("Genome-wide Poisson Threshold (Signal no less than): "+signalBacks.getGenomicModelThreshold("high"));
        if(!noControl)
        	System.out.println("Genome-wide Poisson Threshold (Control no more than): "+ctrlBacks.getGenomicModelThreshold("high"));
        System.out.println("Per-Base Poisson Threshold (Signal): "+signalPerBaseBack.getMaxThreshold('.'));
        if(!noControl)
        	System.out.println("Per-Base Poisson Threshold (Control): "+ctrlPerBaseBack.getMaxThreshold('.'));
        System.out.println("Scaling Factor: "+(noControl ? 1 : control.getScalingFactor()));
        
		//call peak finder for the second time
		testRegions = new ChromosomeGenerator().execute(gen);
		callEnrichedRegions(testRegions, true, false);
		
		// Here we "convert" all peaks to (the more general) features
		signalFeatures.addAll(signalPeaks);
		controlFeatures.addAll(controlPeaks);
		System.out.println("Second round: " + signalPeaks.size()+" enriched signal peaks, "+controlPeaks.size()+" enriched control peaks.\nPrinting peaks to "+getOutName());
		
		if(showGeneAnnotations)
		{
			//Add closest genes
			System.out.println("Adding gene annotations");
			addClosestGenes(signalFeatures);
		}
		
		if(showOtherAnnotations)
		{
			//Add other annotations
			System.out.println("Adding other annotations");
			addRegionAnnotations(signalFeatures);			
		}		
		return signalFeatures;
	}
	
	public void printError(){
		System.err.println("ChipSeqPeakFinder \n");
		this.printArgs();
	}
}
