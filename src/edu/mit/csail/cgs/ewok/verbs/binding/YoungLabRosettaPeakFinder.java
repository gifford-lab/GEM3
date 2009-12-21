/*
 * Created on Mar 20, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.ProbedBindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.MSPProbe;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * @author tdanford
 */
public class YoungLabRosettaPeakFinder implements PeakFinder<MSPProbe> {
        
    private Genome genome;
    private String key;
    private MSPProbe[] triple;
    private double p3Cutoff, pCutoff, nCutoff, twoCutoff;
    
    public YoungLabRosettaPeakFinder(Genome g, String k, double p3, double p1, double n, double two) { 
        genome = g;
        key = k;
        triple = new MSPProbe[3];
        for(int i = 0; i < 3; i++) { triple[i] = null; }
        
        p3Cutoff = p3;
        pCutoff = p1;
        nCutoff = n;
        twoCutoff = two;        
    }

    public YoungLabRosettaPeakFinder(Genome g, String k) {
        genome = g;
        key = k;
        triple = new MSPProbe[3];
        for(int i = 0; i < 3; i++) { triple[i] = null; }
        
        p3Cutoff = 0.001;
        pCutoff = 0.001;
        nCutoff = 0.1;
        twoCutoff = 0.005;
    }
    
    public void reset() { 
        for(int i = 0; i < 3; i++) { triple[i] = null; }        
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Distiller#execute(null)
     */
    public BindingEvent execute(MSPProbe p) {
        for(int i = 0; i < 2; i++) { triple[i] = triple[i+1]; }
        triple[2] = p;
        return distillCurrent();
    }
    
    private boolean currentMeetsPredicate() {
        if(triple[0] == null) { return false; }

		if(!triple[0].containsKey(key)) { 
			//return false; 
			throw new IllegalStateException();
		}
        
        for(int i = 1; i < 3; i++) { 
			if(!triple[i].containsKey(key)) { 
				throw new IllegalStateException();
				//return false; 
			}

            if(!triple[i].getChrom().equals(triple[0].getChrom())) { 
                //return false; 
				throw new IllegalStateException();
            }

			if(triple[i].getLocation() < triple[i-1].getLocation()) { 
				throw new IllegalStateException();
			}
        }

        MSPProbe left = triple[0], middle = triple[1], right = triple[2];
        
        boolean firstCond = triple[2].getPval3(key) <= p3Cutoff;
        if(!firstCond) { return false; }
        
        boolean sc1 = middle.getPval(key) <= pCutoff;
        boolean sc2 = left.getPval(key) <= nCutoff || right.getPval(key) <= nCutoff;
        boolean secondCond = sc1 && sc2;
        if(secondCond) { return true; }
        
        int count = 0;
        for(int i = 0; i < 3; i++) { 
            if(triple[i].getPval(key) <= twoCutoff) { 
                count += 1; 
            }
        }
        boolean thirdCond = count >= 2;
        return thirdCond;
    }
    
    private BindingEvent distillCurrent() { 
        if(currentMeetsPredicate()) { 
            int ns = triple[0].getLocation(), ne = triple[2].getLocation();
            String c = triple[1].getChrom();
            LinkedList<Probe> probes = new LinkedList<Probe>();
            probes.addLast(triple[1]);
            
            return new ProbedBindingEvent(genome, c, ns, ne, triple[1].getRatio(key), triple[1].getPval3(key), key, probes);
        }
        
        return null;        
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Distiller#getCurrent()
     */
    public BindingEvent getCurrent() {
        BindingEvent ret = distillCurrent();
        for(int i = 0; i < 3; i++) { triple[i] = null; }
        return ret;
    }

}
