/*
 * Created on Feb 7, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;
import java.io.*;
import java.text.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.numeric.Numerical;
import edu.mit.csail.cgs.utils.probability.Hypergeometric;

/**
 * @author tdanford
 */
public class BoundaryPValues {
    
    public static class Args { 
        private int N, pos, E;
        
        public Args(int n, int p, int e) { 
            N = n; 
            pos = p;
            E = e;
        }
        
        public int hashCode() { 
            int code = 17;
            code += N; code *= 37;
            code += pos; code *= 37;
            code += E; code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof Args)) { return false; }
            Args a = (Args)o;
            if(N != a.N) { return false; }
            if(pos != a.pos) { return false; }
            if(E != a.E) { return false; }
            return true;
        }
        
        public int compareTo(Args a) { 
            if(N < a.N) { return -1; }
            if(N > a.N) { return 1; }
            if(pos < a.pos) { return -1; }
            if(pos > a.pos) { return 1; }
            if(E < a.E) { return -1; }
            if(E > a.E) { return 1; }
            return 0;
        }
        
        public String toString() { return N + "," + pos + "," + E; }
    }
    
    public static void main(String[] args) { 
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        BoundaryPValues bpv = new BoundaryPValues();
        NumberFormat nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(6);
        
        try { 
            System.out.print(">"); System.out.flush();
            while((line = br.readLine()) != null && line.length() > 0) { 
                StringTokenizer st = new StringTokenizer(line);
                int n = Integer.parseInt(st.nextToken());
                int pos = Integer.parseInt(st.nextToken());
                int e = Integer.parseInt(st.nextToken());
                
                double logpvalue = bpv.getLogPValue(n, pos, e);
                
                System.out.println(n + "," + pos + ":" + e + " --> " + 
                        nf.format(logpvalue) + 
                        " (" + nf.format(Math.exp(logpvalue)) + ")");
                System.out.print(">"); System.out.flush();
            }
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
    }
    
    private Map<Args,Double> pValueCache;
    private ConstrainedChooser chooser;
    private Hypergeometric hypgeom;

    public BoundaryPValues() {
        chooser = new CachingChooser(new GraphChooser());
        hypgeom = new Hypergeometric();
        pValueCache = new HashMap<Args,Double>();
    }
    
    public int getCacheSize() { return pValueCache.size(); }
    public void clearCache() { pValueCache.clear(); }

    public double getLogPValue(int N, int pos, int bestE) {
        if(pos > N) { throw new IllegalArgumentException(); }
        if(bestE > pos || bestE > N-pos) { throw new IllegalArgumentException(bestE + "," + N + "/" + pos); }
        Args a = new Args(N, pos, bestE);
        if(pValueCache.containsKey(a)) { return pValueCache.get(a); }
        
        //double logdenom = hypgeom.log_factorial(N);
        //double logdenom = 0.0;
        double logdenom = hypgeom.log_choose(N, pos);

        Double logsum = null;
        
        int minError = bestE;
        
        for(int e = 0; e <= minError; e++) { 
            if(logsum == null) { 
                logsum = getLogPValueError(N, pos, e, logdenom);
            } else { 
                logsum = Numerical.log_add(logsum, getLogPValueError(N, pos, e, logdenom));
            }
        }

        //System.out.println("Denom: " + logdenom);
        //System.out.println("PValue(" + N + "," + pos + ": " + bestE + ") --> " + Math.exp(logsum));
        
        pValueCache.put(a, logsum);
        return logsum;
    }
    
    private double getLogPValueError(int N, int pos, int error, double denom) {
        Double logsum = null;
        for(int b = 0; b <= N; b++) { 
            if(logsum == null) { 
                logsum = getLogPValueBoundary(N, pos, error, b, denom);
            } else { 
                logsum = Numerical.log_add(logsum, getLogPValueBoundary(N, pos, error, b, denom));
            }
        }
        //System.out.println("\tError(" + error + ") --> " + Math.exp(logsum));
        return logsum;
    }
    
    private double getLogPValueBoundary(int N, int pos, int error, int b, double denom) {
        double logsum = 0.0;
        logsum = Numerical.log_add(getLogPValueDirection(N, pos, error, b, -1, denom),
                    getLogPValueDirection(N, pos, error, b, 1, denom));
        //System.out.println("\t\tBoundary(" + b + ") --> " + Math.exp(logsum));
        return logsum;
    }
    
    private double getLogPValueDirection(int N, int pos, int error, int b, int dir, double denom) {
        if(dir != 1 && dir != -1) { throw new IllegalArgumentException(String.valueOf(dir)); }
        Double logsum = null;
        
        //System.out.println("N: " + N + ", P: " + pos + ", E: " + error + ", b: " + b + ", (" + dir + ")");
        
        if(dir == -1) {
            int numerTest = error - pos + b;
            int param = numerTest / 2;
            if(numerTest % 2 == 0 && param >= 0 && param <= error) { 
                logsum = getLogPValueLeft(N, pos, error, b, numerTest / 2, denom);
            }
        } else { 
            int numerTest = error - pos + (N - b);
            int param = numerTest / 2;
            if(numerTest % 2 == 0 && param >= 0 && param <= error) { 
                logsum = getLogPValueRight(N, pos, error, b, numerTest / 2, denom);
            }
        }
        
        if(logsum != null) { 
            //System.out.println("\t\t\tDirection(" + dir + ") --> " + Math.exp(logsum));
            return logsum;
        } else { 
            //System.out.println("\t\t\tDirection(" + dir + ") --> " + 0.0);
            return -Double.MAX_VALUE;
        }

    }
    
    private double getLogPValueLeft(int N, int pos, int error, int b, int p0, double denom) {
        int p1 = error-p0;
        int s0 = b-p0;
        int s1 = N-b-p1;
        
        //System.out.println("p0/s0: " + p0 + "," + s0 + "  p1/s1: " + p1 + "," + s1);
        
        double left_choose = logOptChoose(b, p0, p1, s1);
        //System.out.println("Left: " + Math.exp(left_choose));
        
        double right_choose = logPOptChoose(N-b, p1, p0, s0);
        //System.out.println("Right: " + Math.exp(right_choose));

        double added = left_choose + right_choose;
        added -= denom;
        
        //System.out.println("Total: " + Math.exp(added));
        return added;        
    }
    
    private double getLogPValueRight(int N, int pos, int error, int b, int p1, double denom) {
        int p0 = error-p1;
        int s1 = N-b-p1;
        int s0 = b-p0;

        //System.out.println("p0/s0: " + p0 + "," + s0 + "  p1/s1: " + p1 + "," + s1);

        double left_choose = logOptChoose(b, p0, p1, s1);        
        //System.out.println("Left: " + Math.exp(left_choose));
        
        double right_choose = logPOptChoose(N-b, p1, p0, s0);
        //System.out.println("Right: " + Math.exp(right_choose));

        double added = left_choose + right_choose;
        added -= denom;
        
        //System.out.println("Total: " + Math.exp(added));
        return added;
    }
    
    private double logOptChoose(int N, int E, int p0, int s0) {
        return chooser.logConstrainedChoose(N, E, p0, s0, true);
    }
    
    private double logPOptChoose(int N, int E, int p0, int s0) { 
        return chooser.logConstrainedChoose(N, E, p0, s0, false);
    }    
    
    public PValuator getPValuator() { return new BoundaryPValuator(); }

    public class BoundaryPValuator implements PValuator {
        
        public BoundaryPValuator() { 
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.PValuator#getLogPValue(edu.mit.csail.psrg.tdanford.boundary.pvalues.BoundaryDataset)
         */
        public PValueResult logPValue(BoundaryDataset ds) {
            int N = ds.size();
            int pos = ds.getNumPositive();
            int errors = ds.getError();
            
            double logPValue = getLogPValue(N, pos, errors);
            PValueResult pvr = 
                new PValueResult(ds.toString(), ds.getError(), ds.getBoundaryDirection(), logPValue);
            return pvr;            
        } 
        
    }
    
}
