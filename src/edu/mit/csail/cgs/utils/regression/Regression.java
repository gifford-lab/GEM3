package edu.mit.csail.cgs.utils.regression;

import cern.colt.matrix.*;
import cern.colt.matrix.linalg.*;
import edu.mit.csail.cgs.utils.Pair;

import java.util.*;

/**
 * <code>Regression</code> provides linear regression (both weighted and unweighted)
 * for a variety of input matrix formats.  If you're looking for linear(double[][], double[]), 
 * have a look at <code>ReuseRegression</code> (which offers some other benefits too)
 *
 */
public class Regression {

    private static final double noiseIncrement = .000001;


    /** The linear() methods that take retrySingular will help mask certain errors.
     *If the indep matrix is singular, the original regression will fail, so linear()
     *  adds a little bit of random noise to each value and tries again.  The idea is that hte noise
     *  will make the matrix non-singular but not change the output too much
     */
    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep, boolean retrySingular) {
        double noisefactor = 0;
        while (noisefactor < 10 * noiseIncrement) {
            try {
                return linear(indep,dep);
            } catch (RuntimeException e) {
                if (!retrySingular) {
                    throw e;
                }
                noisefactor += noiseIncrement;
                for (int i = 0; i < indep.rows(); i++) {
                    for (int j = 0; j < indep.columns(); j++) {
                        indep.setQuick(i,j, indep.getQuick(i,j) + (Math.random() - .5) * noisefactor);
                    }
                }
            }
        }
        throw new IllegalArgumentException("Couldn't do linear regression");
    }
    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep, DoubleMatrix2D weights, boolean retrySingular) {
        double noisefactor = 0;
        while (noisefactor < 10 * noiseIncrement) {
            try {
                return linear(indep,dep,weights);
            } catch (RuntimeException e) {
                if (!retrySingular) {
                    throw e;
                }
                noisefactor += noiseIncrement;
                for (int i = 0; i < indep.rows(); i++) {
                    for (int j = 0; j < indep.columns(); j++) {
                        indep.setQuick(i,j, indep.getQuick(i,j) + (Math.random() - .5) * noisefactor);
                    }
                }
            }
        }
        throw new IllegalArgumentException("Couldn't do linear regression");
    }

    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep) {
        Algebra alg = new Algebra();
        DoubleMatrix2D indeptranspose =alg.transpose(indep); 
         DoubleMatrix2D temp = alg.mult(indeptranspose,indep);
         temp = alg.inverse(temp);
         temp = alg.mult(temp,indeptranspose);
         return alg.mult(temp,dep);
    }

    /* if y is nx1, then weights is nxn.  All elements are zero except weights(i,i) = w_i */
    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep, DoubleMatrix2D weights) {        
        Algebra alg = new Algebra();
        DoubleMatrix2D indeptranspose =alg.transpose(indep); 
        DoubleMatrix2D temp = alg.mult(weights,indep);
        temp = alg.mult(indeptranspose,temp);
        temp = alg.inverse(temp);
        temp = alg.mult(temp,indeptranspose);
        temp = alg.mult(temp,weights);
        return alg.mult(temp,dep);
    }
    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep, DoubleMatrix1D weights) {
        Algebra alg = new Algebra();
        DoubleMatrix2D indeptranspose =alg.transpose(indep); 
        DoubleMatrix2D temp = DoubleFactory2D.dense.make(indep.rows(), indep.columns());
//         System.err.println("indep = " + indep);
//         System.err.println("dep = " + dep);
//         System.err.println("weights= " + weights);
        DoubleMatrix2D newdep = DoubleFactory2D.dense.make(dep.rows(), dep.columns());
        for (int i = 0; i < indep.rows(); i++) {
            double w = weights.getQuick(i);
            newdep.setQuick(i,0, dep.getQuick(i,0) * w);
            for (int j = 0; j < indep.columns(); j++) {
                temp.setQuick(i,j,indep.getQuick(i,j) * w);
            }
        }        
        //        System.err.println("newdep = " + newdep);
        temp = alg.mult(indeptranspose,temp);
        temp = alg.inverse(temp);
        temp = alg.mult(temp,indeptranspose);
        temp = alg.mult(temp,newdep);
        return temp;
    }

    public static DoubleMatrix2D linear(DoubleMatrix2D indep, DoubleMatrix2D dep, double weights[]) {
        Algebra alg = new Algebra();
        DoubleMatrix2D indeptranspose =alg.transpose(indep); 
        DoubleMatrix2D temp = DoubleFactory2D.dense.make(indep.rows(), indep.columns());
        DoubleMatrix2D newdep = DoubleFactory2D.dense.make(dep.rows(), dep.columns());
        for (int i = 0; i < indep.rows(); i++) {
            double w = weights[i];
            newdep.setQuick(i,0, dep.getQuick(i,0) * w);
            for (int j = 0; j < indep.columns(); j++) {
                temp.setQuick(i,j,indep.getQuick(i,j) * w);
            }
        }        
        //        System.err.println("newdep = " + newdep);
        temp = alg.mult(indeptranspose,temp);
        temp = alg.inverse(temp);
        temp = alg.mult(temp,indeptranspose);
        temp = alg.mult(temp,newdep);
        return temp;
    }

    /**
     * Non-Negative Least Squares fitting.  Taken from chapter 23 of Solving Least Squares Problems by Lawson and Hanson 
     * This method uses linear() as a primative.
     */
    public static DoubleMatrix2D nnls(DoubleMatrix2D indep, DoubleMatrix2D dep) {
        boolean verbose = false;

        Algebra alg = new Algebra();
        Set<Integer> P = new TreeSet<Integer>();
        Set<Integer> Z = new TreeSet<Integer>();
        // 1
        DoubleMatrix2D x =  DoubleFactory2D.dense.make(indep.columns(),1);
        for (int i = 0; i < indep.columns(); i++) {
            Z.add(i);
            x.setQuick(i,0,0.0);
        }
        // 2        
        DoubleMatrix2D temp = DoubleFactory2D.dense.make(dep.toArray());
        temp.assign(alg.mult(indep,x), cern.jet.math.Functions.minus);
        DoubleMatrix2D w = alg.mult(alg.transpose(indep),temp);
        // 3
        boolean skiptosix = false;
        if (verbose) {
            System.err.println("trying nnls on ") ;
            for (int i = 0; i < indep.rows(); i++) {
                System.err.print(" " + dep.get(i,0) + " = " );
                for (int j = 0; j < indep.columns(); j++) {
                    System.err.print("  " + indep.get(i,j));
                }
                System.err.println();
            }
        }
        int reps = 0;
        while (Z.size() != 0) {
            if (reps++ > 20 * indep.columns()) {
                throw new IllegalArgumentException("nnls is looping");
            }
            if (verbose) {
                System.err.println("\nP is " + P + "  Z is " + Z + " skiptosix is " + skiptosix + " w=" + w);
                System.err.println(" x is " + x + "\n");
            }
            if (!skiptosix) {
                boolean allnp = true;
                for (int i : Z) {
                    allnp = allnp && w.getQuick(i,0) <= 0;
                }
                if (allnp) { break;}
                double maxval = Double.NEGATIVE_INFINITY;
                // 4
                int t = -1;
                for (int i : Z) {
                    if (w.getQuick(i,0) > maxval) {
                        t = i;
                        maxval= w.getQuick(i,0);
                    }
                }
                //5
                if (t == -1) {
                    throw new RuntimeException("t=-1");
                }
                Z.remove(t);
                P.add(t);
                if (verbose) {
                    System.err.println("Moving " + t + " to P");
                }
            }
            // 6
            DoubleMatrix2D Ep = DoubleFactory2D.dense.make(indep.rows(), P.size(), 0.0);
            DoubleMatrix2D z = null;
            double noiselevel = 0;
            while (z == null) {  
                int cnum = 0;
                for (int j : P) {
                    for (int i = 0; i < Ep.rows(); i++) {
                        Ep.setQuick(i,cnum,indep.getQuick(i,j) + Math.random() * noiselevel);
                    }
                    cnum++;
                }
                try {
                    if (verbose) {
                        System.err.println("trying linreg on ") ;
                        for (int i = 0; i < Ep.rows(); i++) {
                            System.err.print(" " + dep.get(i,0) + " = " );
                            for (int j = 0; j < Ep.columns(); j++) {
                                System.err.print("  " + Ep.get(i,j));
                            }
                            System.err.println();
                        }
                    }

                    z = linear(Ep, dep);
                } catch (IllegalArgumentException e) {
                    noiselevel += .0000001;
                } catch (Exception e) {
                    e.printStackTrace();
                }
                ////System.err.println("Noise level is " + noiselevel);
                if (noiselevel > .00001) {
                    throw new IllegalArgumentException("Couldn't do linreg");
                }

            }
            if (verbose) {
                System.err.println("Solved linreg as " +z);
            }
//             }
            // 7
            boolean allpos = true;
            for (int j = 0; j < P.size(); j++) {
                allpos = allpos && z.getQuick(j,0) > 0;
            }
            if (allpos) {
                int cnum = 0;
                for (int j : Z) {
                    x.setQuick(j,0,0);
                }
                for (int j : P) {
                    x.setQuick(j,0,z.getQuick(cnum,0));
                    cnum++;
                }
                temp = DoubleFactory2D.dense.make(dep.toArray());
                temp.assign(alg.mult(indep,x), cern.jet.math.Functions.minus);
                w = alg.mult(alg.transpose(indep),temp);
                skiptosix = false;
                continue;
            }
            // 8
            double alpha = Double.POSITIVE_INFINITY;
            int q = -1;
            int cnum = 0;
            for (int j : P) {
                if (z.getQuick(cnum,0) > 0) {
                    cnum++;
                    continue;
                }
                double val = x.getQuick(j,0) / (x.getQuick(j,0) - z.getQuick(cnum,0));
                if (val < alpha) {
                    q = j;
                    alpha = val;
                }
                cnum++;
            }
            if (q == -1) {
                throw new RuntimeException("q==-1");
            }
            
            // 9  alpha = minval
            if (verbose) {
                System.err.println("  q=" + q + "  alpha=" + alpha);
            }
            // 10
            for (int j : Z) {
                x.setQuick(j,0, x.getQuick(j,0) * (1 - alpha));
            }
            cnum = 0;
            for (int j : P) {
                x.setQuick(j,0, x.getQuick(j,0) * (1 - alpha) + alpha * z.getQuick(cnum++,0));
            }
            x.setQuick(q,0,0.0);
            if (verbose) {
                System.err.println("  updated x to be : " + x);
            }
            // 11
            ArrayList<Integer> toMove = new ArrayList<Integer>();
            for (int j : P) {
                if (x.getQuick(j,0) == 0) {
                    toMove.add(j);
                }
            }
            if (verbose) {
                System.err.println("Moving " + toMove + " from P to Z");
            }
            P.removeAll(toMove);
            Z.addAll(toMove);
            skiptosix = true;
        }
        if (verbose) {
            System.err.println("Returning " + x);
        }
        return x;
    }

    public static DoubleMatrix2D predict(DoubleMatrix2D indep, DoubleMatrix2D coeffs) {
        Algebra alg = new Algebra();
        return alg.mult(indep,coeffs);
    }

    /**
     * Scores the fit of predicted values to observed values.
     * The returned pair has the 
     *   - score: sum of (actual - obs)^2
     *   - r2
     * where r2 = 1 - ess/tss
     *       ess = sum (actual_i - obs_i)^2
     *       tss = sum (actual_i - mean)^2
     */
    public static Pair<Double,Double> score(double[] predicted, double[] observations) {
        if (predicted.length != observations.length) {
            throw new IllegalArgumentException("Size Mismatch " + predicted.length + 
                                               " vs " + observations.length);
        }
        double mean = 0;
        for (int i = 0; i < observations.length; i++) {
            mean += observations[i];
        }
        mean /= observations.length;
        double ess = 0, tss = 0;
        for (int i = 0; i < predicted.length; i++) {
            ess += Math.pow(observations[i] - predicted[i],2.0);
            tss += Math.pow(observations[i] - mean,2.0);
        }        
        return new Pair(new Double(ess), new Double(1 - ess/tss));
    }
    public static Pair<Double,Double> score(List<Double> predicted, List<Double> observations) {
        if (predicted.size() != observations.size()) {
            throw new IllegalArgumentException("Size Mismatch " + predicted.size() + 
                                               " vs " + observations.size());
        }
        double mean = 0;
        double predmean = 0;
        for (double d : observations) {
            mean += d;
        }
        for (double d: predicted) {
            predmean += d;
        }
        predmean /= predicted.size();        
        mean /= observations.size();
        double ess = 0, tss = 0;
        for (int i = 0; i < predicted.size(); i++) {
            ess += Math.pow(observations.get(i) - predicted.get(i),2.0);
            tss += Math.pow(observations.get(i) - mean,2.0);
        }        
        System.err.println(String.format("mean = %.2f predmean = %.2f, ess = %.2f, tss = %.2f", mean, predmean, ess, tss));
        return new Pair(new Double(ess), new Double(1 - ess/tss));
    }
    public static Pair<Double,Double> score(DoubleMatrix2D predicted, DoubleMatrix2D observations, DoubleMatrix2D weights) {
        Algebra alg = new Algebra();
        if (predicted.rows() != observations.rows()) {
            throw new IllegalArgumentException("Size Mismatch " + predicted.rows() + 
                                               " vs " + observations.rows());
        }
        if (predicted.columns() != 1) {
            throw new IllegalArgumentException("predicted columns != 1");
        }
        if (observations.columns() != 1) {
            throw new IllegalArgumentException("observed columns != 1");
        }

        double mean = 0;
        for (int i = 0; i < observations.rows(); i++) {
            mean += weights.get(i,i) * observations.get(i,0);
        }
        mean /= observations.rows();
        double ess = 0, tss = 0;
        for (int i = 0; i < predicted.rows(); i++) {
            ess += Math.pow(weights.get(i,i) * (observations.get(i,0) - predicted.get(i,0)),2.0);
            tss += Math.pow(weights.get(i,i) * (observations.get(i,0) - mean),2.0);
        }        
        return new Pair(new Double(ess), new Double(1 - ess/tss));
    }

    public static Pair<Double,Double> score(DoubleMatrix2D predicted, DoubleMatrix2D observations, DoubleMatrix1D weights) {
        Algebra alg = new Algebra();
        if (predicted.rows() != observations.rows()) {
            throw new IllegalArgumentException("Size Mismatch " + predicted.rows() + 
                                               " vs " + observations.rows());
        }
        if (predicted.columns() != 1) {
            throw new IllegalArgumentException("predicted columns != 1");
        }
        if (observations.columns() != 1) {
            throw new IllegalArgumentException("observed columns != 1");
        }

        double mean = 0;
        for (int i = 0; i < observations.rows(); i++) {
            mean += weights.get(i) * observations.get(i,0);
        }
        mean /= observations.rows();
        double ess = 0, tss = 0;
        for (int i = 0; i < predicted.rows(); i++) {
            ess += Math.pow(weights.get(i) * (observations.get(i,0) - predicted.get(i,0)),2.0);
            tss += Math.pow(weights.get(i) * (observations.get(i,0) - mean),2.0);
        }        
        return new Pair(new Double(ess), new Double(1 - ess/tss));
    }

    public static Pair<Double,Double> score(DoubleMatrix2D predicted, DoubleMatrix2D observations) {
        Algebra alg = new Algebra();
        if (predicted.rows() != observations.rows()) {
            throw new IllegalArgumentException("Size Mismatch " + predicted.rows() + 
                                               " vs " + observations.rows());
        }
        if (predicted.columns() != 1) {
            throw new IllegalArgumentException("predicted columns != 1");
        }
        if (observations.columns() != 1) {
            throw new IllegalArgumentException("observed columns != 1");
        }

        double mean = 0;
        for (int i = 0; i < observations.rows(); i++) {
            mean += observations.get(i,0);
        }
        mean /= observations.rows();
        double ess = 0, tss = 0;
        for (int i = 0; i < predicted.rows(); i++) {
            ess += Math.pow(observations.get(i,0) - predicted.get(i,0),2.0);
            tss += Math.pow(observations.get(i,0) - mean,2.0);
        }        
        return new Pair(new Double(ess), new Double(1 - ess/tss));
    }

}

