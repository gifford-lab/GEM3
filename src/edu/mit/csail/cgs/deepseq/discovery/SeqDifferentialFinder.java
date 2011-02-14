package edu.mit.csail.cgs.deepseq.discovery;

import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.deepseq.*;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class SeqDifferentialFinder extends DifferentialEnrichmentFinder{

	public SeqDifferentialFinder(DeepSeqExpt signal) {
		super(signal, null);
	}
	public SeqDifferentialFinder(DeepSeqExpt signal, DeepSeqExpt ctrl) {
		super(signal, ctrl);
	}
	public SeqDifferentialFinder(String[] args){
		super(args);
	}
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		
		System.out.println("Welcome to the SeqDifferentialFinder!!\n\n");
		
		long start = System.currentTimeMillis();
		SeqDifferentialFinder finder = new SeqDifferentialFinder(args);
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
	
	
	public List<Feature> execute(){
		//Initialize background models
		initializeBackgrounds();
		
		//Print some more info
		System.out.println("Sliding Window: Size="+binWidth+" Offset="+binStep);
        System.out.println("Read Length: "+readLength);
        System.out.println("Genome-wide Poisson Threshold (Signal no less than): "+signalBacks.getGenomicModelThreshold("high"));
        if(!noControl)
        	System.out.println("Genome-wide Poisson Threshold (Control no more than): "+ctrlBacks.getGenomicModelThreshold("high"));
//        System.out.println("Per-Base Poisson Threshold (Signal): "+signalPerBaseBack.getMaxThreshold('.'));
//      if(!noControl)
//        	System.out.println("Per-Base Poisson Threshold (Control): "+ctrlPerBaseBack.getMaxThreshold('.'));
        System.out.println("Scaling Factor: "+(noControl ? 1 : control.getScalingFactor()));
        
		//Initialize the iterator for the test regions/genes
		Iterator<Region> testRegions = new ChromosomeGenerator().execute(gen);
		
		//Call the peak finder for the first round (finds enriched regions for the meta-peak and sets the control scaling factor
		callEnrichedRegions(testRegions, true, false);
		System.out.println("First round: "+signalEvents.size()+" enriched signal peaks, "+controlEvents.size()+" enriched control peaks.");
				
		// Here we "convert" all peaks to (the more general) features
		signalFeatures.addAll(signalEvents);
		controlFeatures.addAll(controlEvents);
		
		if(showGeneAnnotations)
		{	//Add closest genes
			System.out.println("Adding gene annotations");
			addClosestGenes(signalFeatures);
		}
		
		if(showOtherAnnotations)
		{	//Add other annotations
			System.out.println("Adding other annotations");
			addRegionAnnotations(signalFeatures);			
		}		
		return signalFeatures;
	}
	
	public void printError(){
		System.err.println("SeqDifferentialFinder \n");
		this.printArgs();
	}
}
