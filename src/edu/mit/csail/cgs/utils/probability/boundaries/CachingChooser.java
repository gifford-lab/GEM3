/*
 * Created on Feb 10, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class CachingChooser implements ConstrainedChooser {
    
    private Map<Arguments,Double> optMap, pOptMap;
    private ConstrainedChooser chooser;
    
    public CachingChooser(ConstrainedChooser c) {
        optMap = new TreeMap<Arguments,Double>();
        pOptMap = new TreeMap<Arguments,Double>();
        chooser = c;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.ConstrainedChooser#logConstrainedChoose(int, int, int, int, boolean)
     */
    public double logConstrainedChoose(int N, int E, int p0, int s0, boolean optimal) {
        
        Arguments a = new Arguments(N, E, p0, s0);
        
        if(!optMap.containsKey(a)) { 
            optMap.put(a, chooser.logConstrainedChoose(N, E, p0, s0, true));
            pOptMap.put(a, chooser.logConstrainedChoose(N, E, p0, s0, false));
        }
        
        if(optimal) { 
            return optMap.get(a);
        } else { 
            return pOptMap.get(a);
        }
    }

    private static class Arguments implements Comparable<Arguments> {
        
        public int N, E, p0, s0;
        
        public Arguments(int N, int E, int p0, int s0) { 
            this.N = N; this.E = E;
            this.p0 = p0; this.s0 = s0;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Arguments)) { return false; }
            Arguments a = (Arguments)o;
            if(N != a.N) { return false; }
            if(E != a.E) { return false; }
            if(p0 != a.p0) { return false; }
            if(s0 != a.s0) { return false; }
            return true;
        }
        
        public int hashCode() { 
            int code = 17;
            code += N; code *= 37; 
            code += E; code *= 37;
            code += p0; code *= 37;
            code += s0; code *= 37;
            return code;
        }

        /* (non-Javadoc)
         * @see java.lang.Comparable#compareTo(java.lang.Object)
         */
        public int compareTo(Arguments a) {
            if(N < a.N) { return -1; }
            if(N > a.N) { return 1; }
            if(E < a.E) { return -1; }
            if(E > a.E) { return 1; }
            if(p0 < a.p0) { return -1; }
            if(p0 > a.p0) { return 1; }
            if(s0 < a.s0) { return -1; }
            if(s0 > a.s0) { return 1; }
            return 0;
        }
    }
}
