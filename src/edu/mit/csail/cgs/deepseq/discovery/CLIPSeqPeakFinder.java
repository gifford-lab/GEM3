package edu.mit.csail.cgs.deepseq.discovery;

import java.sql.SQLException;
import java.util.List;

import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class CLIPSeqPeakFinder extends StatisticalPeakFinder{
	
	public CLIPSeqPeakFinder(DeepSeqExpt signal) {
		super(signal, null);
		read5PrimeExt = 0;
		read3PrimeExt=0;
		readShift=0;
		signal.setFivePrimeExt((int)read5PrimeExt); 
		signal.setThreePrimeExt((int)read3PrimeExt); 
		signal.setShift((int)readShift); 
		if(!noControl){
			control.setFivePrimeExt((int)read5PrimeExt);
			control.setThreePrimeExt((int)read3PrimeExt);
			control.setShift((int)readShift);
		}
	}
	public CLIPSeqPeakFinder(DeepSeqExpt signal, DeepSeqExpt ctrl) {
		super(signal, ctrl);
		read5PrimeExt = 0;
		read3PrimeExt=0;
		readShift=0;
		signal.setFivePrimeExt((int)read5PrimeExt); 
		signal.setThreePrimeExt((int)read3PrimeExt); 
		signal.setShift((int)readShift); 
		if(!noControl){
			control.setFivePrimeExt((int)read5PrimeExt);
			control.setThreePrimeExt((int)read3PrimeExt);
			control.setShift((int)readShift);
		}
	}
	public CLIPSeqPeakFinder(String[] args){
		super(args);
		//Reset read extensions with new defaults
		read5PrimeExt = Args.parseDouble(args,"read5ext",0);
		read3PrimeExt = Args.parseDouble(args,"read3ext",0);
		readShift = Args.parseDouble(args,"readshift",0);
		signal.setFivePrimeExt((int)read5PrimeExt); 
		signal.setThreePrimeExt((int)read3PrimeExt); 
		signal.setShift((int)readShift); 
		if(!noControl){
			control.setFivePrimeExt((int)read5PrimeExt);
			control.setThreePrimeExt((int)read3PrimeExt);
			control.setShift((int)readShift);
		}
		
		//can also default dbacks here
		dbacks= (List<Integer>) Args.parseIntegers(args, "dynback");
		//if(dbacks.size()==0)
		//	dbacks.add(new Integer(5000));
	}
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		System.out.println("Welcome to the CLIPSeqPeakFinder\nLoading experiments...");
		CLIPSeqPeakFinder finder = new CLIPSeqPeakFinder(args);
		finder.setStrandedFinding(true);
		finder.setNeedleFilter(false);
		finder.setTowerFilter(false);

		System.out.println("Finding enriched regions");
		List<Feature> peaks = finder.execute();
		System.out.println(peaks.size()+" enriched regions found.\nPrinting peaks to "+finder.getOutName()+"_signal.peaks");
		finder.printFeatures();
		System.out.println("Printing sequences to "+finder.getOutName()+"_signal.seq");
		finder.printPeakSequences();
		finder.cleanup();
	}

	public void printError(){
		System.err.println("CLIPSeqPeakFinder \n");
		this.printArgs();
	}
}
