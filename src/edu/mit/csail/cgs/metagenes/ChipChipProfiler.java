package edu.mit.csail.cgs.metagenes;

import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ChipChipProfiler  implements PointProfiler<Point, Profile>{

	private Genome gen;
	private BinningParameters params=null;
	private ChipChipLocator loc=null;
	private boolean includeIP=true, includeWCE=true;
	
	public ChipChipProfiler(BinningParameters bp, Genome g, ChipChipLocator experiment, boolean useIP, boolean useWCE){
		this(bp, g, experiment);
		includeIP=useIP;
		includeWCE=useWCE;
	}
	public ChipChipProfiler(BinningParameters bp, Genome g, ChipChipLocator experiment){
		gen=g;
		params=bp; 
		loc = experiment;
	}
	
	public BinningParameters getBinningParameters() {
		return params;
	}

	public Profile execute(Point a) {
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		Region query = new Region(gen, a.getChrom(), start, end);
	
		double [] interp = new double[window];
		for(int i=0; i<interp.length; i++){interp[i]=0;}
		
		try {
			ChipChipData data = loc.createObject();
			data.window(a.getChrom(), start, end);
		
			double[] datapoints = new double[data.getCount()];
			int[] positions = new int[data.getCount()];
			for(int i = 0; i < data.getCount(); i++) {
				positions[i] = data.getPos(i);
				double val=0;
				if(includeIP && includeWCE)
					val = getMeanRatio(data, i);
				else if(!includeIP && includeWCE)
					val = getMeanWCE(data, i);
				else
					val = getMeanIP(data, i);
				
				if(Double.isInfinite(val) || Double.isNaN(val))
					datapoints[i]=0;
				else
					datapoints[i]=val;
			}
			
			double lastVal=0;
			int lastOff=0;
			for(int i = 0; i < data.getCount(); i++) {
				int offset = Math.max(0, positions[i]-query.getStart());
				double span = (double)(offset-lastOff);
				double inc = (datapoints[i]-lastVal)/span;
				double count=0;
				for(int o=lastOff; o<offset; o++){
					if(o<interp.length)
						interp[o] = lastVal+(count*inc);
					count++;
				}
				lastOff=offset;
				lastVal = datapoints[i];
			}
			
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//Transform back to array of bins
		int b=0;
		for(int i=0; i<interp.length; i+=params.getBinSize()){
			double sum=0, count=0;
			for(int j=i; j<i+params.getBinSize() && j<interp.length; j++){
				if(!strand)
					sum+=interp[window-j-1];
				else
					sum+=interp[j];
				count++;
			}
			if(b<array.length)
				array[b]=sum/count;
			b++;
		}
		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	private static double getMeanRatio(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(!Double.isInfinite(d.getRatio(i, j)) && !Double.isNaN(d.getRatio(i, j))){
				sum += d.getRatio(i, j);
				count += 1;
			}
		}
		return count > 0 ? sum / count : 0.0;
	}
	private static double getMeanIP(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(!Double.isInfinite(d.getIP(i, j)) && !Double.isNaN(d.getIP(i, j))){
				sum += d.getIP(i, j);
				count += 1;
			}
		}
		return count > 0 ? sum / count : 0.0;
	}
	private static double getMeanWCE(ChipChipData d, int i) { 
		double sum = 0.0;
		double count = 0;
		for(int j = 0; j < d.getReplicates(i); j++) { 
			if(!Double.isInfinite(d.getWCE(i, j)) && !Double.isNaN(d.getWCE(i, j))){
				sum += d.getWCE(i, j);
				count += 1;
			}
		}
		return count > 0 ? sum / count : 0.0;
	}
	
	//No cleanup
	public void cleanup(){}
}
