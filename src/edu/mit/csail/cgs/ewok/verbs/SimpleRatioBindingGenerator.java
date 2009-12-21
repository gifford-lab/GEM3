package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Region;

/**
 * @author tdanford
 *
 */
public class SimpleRatioBindingGenerator<RegionType extends Region> implements Expander<RegionType,BindingEvent>{
	
    private ChipChipData data;
    private double thresh;
    private int maxExtend;
    private String type;

    public SimpleRatioBindingGenerator(ChipChipData data, double thresh1) { 
        this.data = data;
        this.thresh = thresh1;
        maxExtend = 1000;
        type = data.getName();
    }
    
    public Iterator<BindingEvent> execute(RegionType r) {
    	LinkedList<BindingEvent> extents = new LinkedList<BindingEvent>();
    	try {
			data.window(r.getChrom(), r.getStart(), r.getEnd());
			
			int start = -1;
			int ppos = -1;
			double maxRatio = 0.0;
			double meanRatio = 0.0;
			int probeCount = 0;
			int pos = -1;
			
			for(int i = 0; i < data.getCount(); i++) { 
				double ratio = getMeanRatio(i);
				pos = data.getPos(i);
				
				if(ratio >= thresh) { 
					if(ppos == -1) {
						maxRatio = ratio;
						meanRatio = ratio;
						probeCount = 1;
						if(i == 0) { 
							ppos = pos - maxExtend;
						} else { 
							int prevPosition = data.getPos(i-1);
							ppos = Math.max(prevPosition, pos-maxExtend);
						}
					} else { 
						maxRatio = Math.max(maxRatio, ratio);
						probeCount += 1;
						meanRatio += ratio;
					}
				} else { 
					if(ppos != -1) {
						meanRatio /= (double)probeCount;
						int end = pos + maxExtend;
						if(i < data.getCount()-1) {  
							end = Math.min(end, data.getPos(i+1));
						}
						BindingEvent ext = new BindingEvent(r.getGenome(), r.getChrom(), ppos, end, maxRatio, meanRatio, type);
						extents.add(ext);
					}
				}
			}
			
			if(ppos != -1) {
				meanRatio /= (double)probeCount;
				int end = pos + maxExtend;
				BindingEvent ext = new BindingEvent(r.getGenome(), r.getChrom(), ppos, end, maxRatio, meanRatio, type);
				extents.add(ext);
			}
			
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		return extents.iterator();
    }
    
    private double getMeanRatio(int idx) { 
    	double sum = 0.0;
    	int count = data.getReplicates(idx);
    	
    	for(int j = 0; j < count; j++) { 
    		sum += data.getRatio(idx, j);
    	}
  
    	return sum / (double)count;
    }

}
