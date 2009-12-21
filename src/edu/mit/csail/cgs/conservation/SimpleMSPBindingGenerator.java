package edu.mit.csail.cgs.conservation;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMSP;
import edu.mit.csail.cgs.datasets.general.Region;

public class SimpleMSPBindingGenerator implements Expander<Region,BindingExtent> {

    private ChipChipMSP data;
    private double pCutoff;
    private int cutoffDistance;

    public SimpleMSPBindingGenerator(ChipChipMSP data) { 
        this.data = data;
        pCutoff = 0.001;
        cutoffDistance = 1000;
    }
    
    public SimpleMSPBindingGenerator(ChipChipMSP data, double p) { 
        this.data = data;
        pCutoff = p;
        cutoffDistance = 1000;
    }
    
    public boolean meetsCriteria(ChipChipMSP data, int i) {
    	return data.getPval(i) <= pCutoff;
    }
    
    public Iterator <BindingExtent> execute (Region r) {
        LinkedList<BindingExtent> results = new LinkedList<BindingExtent>();
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            int start = -1;
            int end = -1;
            int current = -1;
            double maxProb = 0.0, maxRatio = 0.0;

            for (int i = 1; i < data.getCount()-1; i++) {
            	current = data.getPos(i);
            	
            	// this first check is to make sure we haven't traversed a "break" in the
            	// probe tiling, in which case if we're in the "middle" of a peak, we just
            	// tie it off and call it right here.
            	if(start != -1 && current-data.getPos(i-1) > cutoffDistance) {
                    results.add(new BindingExtent(r.getGenome(),
                            r.getChrom(),
                            start, end,
                            maxRatio, maxProb, "MSP", start, end));
                    start = -1;
                    maxProb = maxRatio = 0.0;
            	}

            	if(meetsCriteria(data, i)) { 
                    if(start == -1) { 
                        start = Math.max(data.getPos(i-1), current-cutoffDistance);                        
                        maxProb = data.getPval(i);
                        maxRatio = data.getMedianOfRatios(i);
                    } else { 
                        maxProb = Math.max(maxProb, data.getPval(i));
                        maxRatio = Math.max(maxRatio, data.getMedianOfRatios(i));
                    }
                    end = Math.min(data.getPos(i+1), current+cutoffDistance);
                } else { 
                    if(start != -1) { 
                        results.add(new BindingExtent(r.getGenome(),
                                r.getChrom(),
                                start, end,
                                maxRatio, maxProb, "MSP", start, end));
                        start = -1;
                        maxProb = maxRatio = 0.0;
                    }
                }
            }
            
            if(start != -1) {
                results.add(new BindingExtent(r.getGenome(),
                        r.getChrom(),
                        start, end,
                        maxRatio, maxProb, "MSP", start, end));                
            }
            
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();        
    }
}
