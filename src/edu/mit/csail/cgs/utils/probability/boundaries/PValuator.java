/*
 * Created on Feb 24, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.numeric.Numerical;
import edu.mit.csail.cgs.utils.probability.*;

/**
 * @author tdanford
 */
public interface PValuator {
    
    public PValueResult logPValue(BoundaryDataset ds);
    
    public static class MannWhitney implements PValuator {
        
        private int[] argArray;
        private MannWhitneyEquation eq;
        
        public MannWhitney() { 
            argArray = new int[3];
            eq = new MannWhitneyEquation();
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.PValuator#logPValue(edu.mit.csail.psrg.tdanford.boundary.pvalues.BoundaryDataset)
         */
        public PValueResult logPValue(BoundaryDataset ds) {
            
            int n = 0, m = 0, u = 0;
            int dir = ds.getBoundaryDirection();
            String dsRep = null;
            
            if(dir == 1) { 
                n = ds.getNumNegative();
                m = ds.getNumPositive();
                dsRep = ds.createStringRep(1);
            } else { 
                n = ds.getNumPositive();
                m = ds.getNumNegative();
                dsRep = ds.createStringRep(0);
            }
            //u = MannWhitneyEquation.calculateUFromT(m, n, dsRep);
            u = MannWhitneyEquation.calculateU(dsRep);
            
            
            double pv = eq.getLowerPValue(m, n, u);

            System.out.println("[" + m + "," + n + "," + u + "] --> " + pv);
            
            double logPV = Math.log(pv);
            
            PValueResult pvr =  
                new PValueResult(ds.toString(), ds.getError(), dir, logPV);
            
            return pvr;
        } 
        
    }

    public static class KolmogorovSmirnov implements PValuator {
        
        public KolmogorovSmirnov() { 
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.PValuator#logPValue(edu.mit.csail.psrg.tdanford.boundary.pvalues.BoundaryDataset)
         */
        public PValueResult logPValue(BoundaryDataset ds) {
            int N = ds.size();
            int pos = ds.getNumPositive();
            int errors = ds.getError();
            
            double[] posa = ds.getPositiveArray();
            double[] nega = ds.getNegativeArray();
            Pair<Double,Double> ks = Numerical.kstwo(posa, nega);
            double logPValue = Math.log(ks.getLast());            
            
            PValueResult pvr = 
                new PValueResult(ds.toString(), ds.getError(), ds.getBoundaryDirection(), logPValue);
            return pvr;            

        } 
        
    }
    
    public static class Spearman implements PValuator {
        
        public Spearman() { 
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.PValuator#logPValue(edu.mit.csail.psrg.tdanford.boundary.pvalues.BoundaryDataset)
         */
        public PValueResult logPValue(BoundaryDataset ds) {
            int N = ds.size();
            int pos = ds.getNumPositive();
            int errors = ds.getError();
            
            double[] posa = ds.getPositiveArray();
            double[] nega = ds.getNegativeArray();
            
            Numerical.SpearmanResult res = Numerical.spear(posa, nega);
            double logPValue = Math.log(res.getProbRS());
            int dir = res.getRS() >= 0.0 ? 1 : -1;
            
            PValueResult pvr = 
                new PValueResult(ds.toString(), ds.getError(), dir, logPValue);
            
            return pvr;            

        } 
        
    }
}
