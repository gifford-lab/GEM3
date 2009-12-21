package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMSP;
import edu.mit.csail.cgs.datasets.general.Region;

public class CustomMSPBindingGenerator implements Expander<Region,BindingExtent> {

    private ChipChipMSP data;
    private double p3Cutoff, pCutoff, nCutoff, twoCutoff;
    private int cutoffDistance;

    public CustomMSPBindingGenerator(ChipChipMSP data) { 
        this.data = data;
        p3Cutoff = 0.001;
        pCutoff = 0.001;
        nCutoff = 0.1;
        twoCutoff = 0.005;
        cutoffDistance = 1000;
    }
    
    public CustomMSPBindingGenerator(ChipChipMSP data, double p3, double p, double n, double two) { 
        this.data = data;
        p3Cutoff = p3;
        pCutoff = p;
        nCutoff = n;
        twoCutoff = two;
        cutoffDistance = 1000;
    }
    
    public boolean meetsCriteria(ChipChipMSP data, int i) { 
        boolean firstCond = data.getPval3(i) <= p3Cutoff;
        if(firstCond) { return true; }
        
        boolean sc1 = data.getPval(i) <= pCutoff;
        boolean sc2 = data.getPval(i-1) <= nCutoff || data.getPval(i+1) <= nCutoff;
        boolean secondCond = sc1 && sc2;
        if(secondCond) { return true; }
        
        int count = 0;
        for(int j = -1; j <= 1; j++) { 
            if(data.getPval(i+j) <= twoCutoff) { 
                count += 1; 
            }
        }
        boolean thirdCond = count >= 2;
        return thirdCond;
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
                        maxProb = data.getPval3(i);
                        maxRatio = data.getMedianOfRatios(i);
                    } else { 
                        maxProb = Math.max(maxProb, data.getPval3(i));
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
