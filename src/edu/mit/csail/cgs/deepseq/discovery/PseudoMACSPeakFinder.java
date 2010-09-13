package edu.mit.csail.cgs.deepseq.discovery;

import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.BindingModelGenerator;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;

public class PseudoMACSPeakFinder extends StatisticalPeakFinder{

	public PseudoMACSPeakFinder(DeepSeqExpt signal) {
		super(signal, null);
	}
	public PseudoMACSPeakFinder(DeepSeqExpt signal, DeepSeqExpt ctrl) {
		super(signal, ctrl);
	}
	public PseudoMACSPeakFinder(String[] args){
		super(args);
		
	}
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		System.out.println("Welcome to the pseudoMACSPeakFinder\nLoading experiments...");
		PseudoMACSPeakFinder finder = new PseudoMACSPeakFinder(args);
		finder.setStrandedFinding(false);
		System.out.println("Finding enriched regions");
		List<Feature> peaks = finder.execute(); 
		finder.printFeatures();
		finder.printFeatures(false);
		System.out.println("Output file root: "+finder.getOutName());
		finder.printPeakSequences();
		finder.cleanup();
	}
	
	public List<Feature> execute(){
				
		//Initialize background models
		initializeBackgrounds();
		
		//Print some more info
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
		callEnrichedRegions(testRegions, true, false);
		System.out.println("First round of peak-calling: "+signalPeaks.size()+" signal peaks, "+controlPeaks.size()+" control peaks.");
		
		//make the meta-peak
		if(buildBindingModel){
			BindingModelGenerator bmg = new BindingModelGenerator(gen, signal, signalPeaks, controlPeaks);
			bmg.setPFilter(1e-8);
			bmg.setOverRepFilter(10);
			setBindingModel(bmg.execute()); //Update shift, ext
			bindingModel.smooth(BindingModel.SMOOTHING_STEPSIZE, BindingModel.SMOOTHING_AVG_PTS);
			setModelPeaks(true);
			bmg.printModel(new String(getOutName()+"_model.txt"));
		}
				
		signalFeatures.addAll(signalPeaks);
		controlFeatures.addAll(controlPeaks);
		System.out.println(signalPeaks.size()+" enriched regions found in signal, "+controlPeaks.size()+" in control.\nPrinting peaks to "+getOutName());
		
		//Add closest genes
		System.out.println("Adding gene annotations");
		addClosestGenes(signalFeatures);
		
		//Add other annotations
		System.out.println("Adding other annotations");
		addRegionAnnotations(signalFeatures);
		
		return signalFeatures;
	}
	
	public void printError(){
		System.err.println("pseudoMACSPeakFinder \n" +
				"");
	}
}
