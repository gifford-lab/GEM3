package edu.mit.csail.cgs.utils.regression;

import cern.colt.matrix.*;
import cern.colt.matrix.linalg.*;

/**
 * Does one, two, or three variable linear regression.  Caches the inputs (X, Y, W) and
 *  can perform each regression on only a subset of the data.  This can be much more efficient than
 *  the routines in Regression if you're working on many subsets of the same inputs.
 *
 *  Solving weighted linear regression is
 *  (X' W X)^-1 WY
 *  The regression routines include a noisefactor: if the matrix inversion fails because the matrix is
 *  singular, ReuseRegression adds a small amount of random noise to X' W X and tries again; it does
 *  this ten times before giving up.  This means that some solutions may be approximate, but it avoids the headaches of 
 *  methods failing.  
 *
 *  ReuseRegression will throw an IllegalArgumentException if the inputs contain any NaN values.
 */

public class ReuseRegression {
    private int COLUMNS = 3;
    private static final double noiseIncrement = .00000001;
    private cern.jet.random.Normal normal;
    private double[][] X;
    private double[] weights, WY, XTWY, Y;
    private int n;

    /**
     * Sets up a new linear regression solver with equally weighted points.  This
     * class does not add a constant to X, so you need to include a column
     * of 1s if you want it.
     *
     * @param X is an nxm matrix.
     * @param Y is an nx1 vector
     */
    public ReuseRegression (double[][] X, double[] Y) {
        this.X = X;
        this.Y = Y;
        this.weights = new double[Y.length];
        for (int i =0; i < weights.length; i++) {
            weights[i] = 1;
        }
        init();
    }

    /**
     * Sets up a new linear regression solver. This
     * class does not add a constant to X, so you need to include a column
     * of 1s if you want it.
     *
     * @param X is an nxm matrix.
     * @param Y is an nx1 vector.
     * @param W is an nx1 vector of weights.
     */
    public ReuseRegression (double[][] X, double[] Y, double[] weights) {
        this.X = X;
        this.Y = Y;
        this.weights = weights;
        init();
    }
    private void init() {
        COLUMNS = X[0].length;
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.DRand();
        normal = new cern.jet.random.Normal(0.0,1.0,engine);
        n = X.length;
        WY = new double[n];
        XTWY = new double[COLUMNS];        
        for (int i = 0; i < n; i++) {
            WY[i] = Y[i] * weights[i];            
        }

        for (int i = 0; i < n; i++) {
            if (Double.isNaN(Y[i])) {
                throw new IllegalArgumentException("NaN at index " + i + " of Y");
            }
            for (int j = 0; j < X[i].length; j++) {
                if (Double.isNaN(X[i][j])) {
                    throw new IllegalArgumentException("NaN at index " + i + "," + j + "  of X");
                }   
            }
        }
    }

    /**
     * returns a 1-element array containing c where
     * c minimizes
     *  sum_i=1^n w_i * (y_i - c) ^ 2
     */
    public double[] fitOneVar(int first, int last) {
        double sum = 0, sumweights = 0;        
        for (int i = first; i <= last; i++) {            
            sum += WY[i];
            sumweights += weights[i];
        }
        double[] out = new double[1];
        out[0] = sum / sumweights;
        return out;
    }

    /**
     * returns a 2-element array containing a and b
     * to minimize
     *  sum_i=1^n w_i * (y_i - (x_i0 * a + x_i1 * b) ^ 2
     *
     * The implementation unrolls all the matrix operations
     * and includes a hardcoded formula for 2x2 matrix inversion.
     */
    public double[] fitTwoVar(int first, int last) {
        double tIA, tIB, tIC, tID;
        double a, b, c, d;
        double[] output = new double[2];
        tIA = 0; tIB = 0;
        tIC = 0; tID = 0;
        XTWY[0] = 0; XTWY[1] = 0;
        
        for (int i = first; i <= last; i++) {
            double wy = WY[i];
            double weight = weights[i];

            XTWY[0] += wy * X[i][0];
            XTWY[1] += wy * X[i][1];

            tIA += weight * X[i][0] * X[i][0];
            tIB += weight * X[i][1] * X[i][0];
            tID += weight * (X[i][1] * X[i][1]);
        }
        tIC = tIB;
        double noisefactor = 0;
        while (noisefactor < 10 * noiseIncrement) {
            double det = (tIA * tID - tIC * tIB);
            if (det == 0) {
                //                System.err.println("DET=0 noisefactor="+noisefactor);
                noisefactor += noiseIncrement;
                tIA += tIA * noisefactor * normal.nextDouble();
                tIB += tIB * noisefactor * normal.nextDouble();
                tIC += tIC * noisefactor * normal.nextDouble();
                tID += tID * noisefactor * normal.nextDouble();
                continue;
            }
            if (Double.isNaN(det)) {
                System.err.println("det is NaN from " + tIA + ", " + tIB + ", " + tIC + ", " + tID);
                continue;
            }

            a = tID / det;
            b = -1 * tIB / det;
            c = -1 * tIC / det;
            d = tIA / det;
            if (Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(c) || Double.isNaN(d)) {
                System.err.println("a=" + a + "  b=" + b + "  c=" + c + "  d=" + d);
                continue;
            }

            output[0] = a * XTWY[0] + b * XTWY[1];
            output[1] = c * XTWY[0] + d * XTWY[1];
            if (Double.isNaN(output[0]) || Double.isNaN(output[1])) {
                System.err.println("a=" + a + "  b=" + b + "  c=" + c + "  d=" + d);
                System.err.println("output[0] = " + output[0] + "  output[1]=" + output[1]);
                continue;
            }

            return output;
        }
        throw new IllegalArgumentException("Couldn't do regression");
    }

    /**
     * returns a 3-element array containing a, b, and c
     * to minimize
     *  sum_i=1^n w_i * (y_i - (x_i0 * a + x_i1 * b + x_i2 * c) ^ 2
     *
     * The implementation unrolls all the matrix operations
     * and includes a hardcoded formula for 3x3 matrix inversion.
     */
    public double[] fitThreeVar(int first, int last) {
        /* tI == to Invert
           i = inverted */
        double tI00, tI01, tI02, tI10, tI11, tI12, tI20, tI21, tI22;
        double i00, i01, i02, i10, i11, i12, i20, i21, i22;

        tI00 = 0; tI01 = 0; tI02 = 0; 
        tI10 = 0; tI11 = 0; tI12 = 0; 
        tI20 = 0; tI21 = 0; tI22 = 0; 
        double[] output = new double[3];
        XTWY[0] = 0; XTWY[1] = 0; XTWY[2] = 0;
        for (int i = first; i <= last; i++) {
            double wy = WY[i];
            double weight = weights[i];

            XTWY[0] += wy * X[i][0];
            XTWY[1] += wy * X[i][1];
            XTWY[2] += wy * X[i][2];

            tI00 += weight * X[i][0] * X[i][0];
            tI01 += weight * X[i][0] * X[i][1];
            tI02 += weight * X[i][0] * X[i][2];
            tI12 += weight * X[i][1] * X[i][2];
            tI22 += weight * X[i][2] * X[i][2];
            tI11 += weight * X[i][1] * X[i][1];
        }
        tI10 = tI01;
        tI20 = tI02;
        tI21 = tI12;

        double noisefactor = 0;
        while (noisefactor < 10 * noiseIncrement) {
            double det = tI00 * (tI11 * tI22 - tI12 * tI21) -  tI01 * (tI10 * tI22 - tI12 * tI20) + tI02 * (tI10 * tI21 - tI11 * tI20);
            if (det == 0) {
                //                System.err.println("DET=0, noisefactor=" + noisefactor);
                noisefactor += noiseIncrement;
                tI00 += tI00 * noisefactor * normal.nextDouble();
                tI10 += tI10 * noisefactor * normal.nextDouble();
                tI20 += tI20 * noisefactor * normal.nextDouble();
                tI01 += tI01 * noisefactor * normal.nextDouble();
                tI11 += tI11 * noisefactor * normal.nextDouble();
                tI21 += tI21 * noisefactor * normal.nextDouble();
                tI02 += tI02 * noisefactor * normal.nextDouble();
                tI12 += tI12 * noisefactor * normal.nextDouble();
                tI22 += tI22 * noisefactor * normal.nextDouble();
                continue;
            }

            det = 1.0 / det;
            i00 = det * (tI11 * tI22 - tI12 * tI21);
            i01 = det * (tI02 * tI21 - tI01 * tI22);
            i02 = det * (tI01 * tI12 - tI02 * tI11);

            i10 = det * (tI12 * tI20 - tI10 * tI22);
            i11 = det * (tI00 * tI22 - tI02 * tI20);
            i12 = det * (tI02 * tI10 - tI00 * tI12);

            i20 = det * (tI10 * tI21 - tI11 * tI20);
            i21 = det * (tI01 * tI20 - tI00 * tI21);
            i22 = det * (tI00 * tI11 - tI01 * tI10);
            
            output[0] = i00 * XTWY[0] + i01 * XTWY[1] + i02 * XTWY[2];
            output[1] = i10 * XTWY[0] + i11 * XTWY[1] + i12 * XTWY[2];
            output[2] = i20 * XTWY[0] + i21 * XTWY[1] + i22 * XTWY[2];
            return output;

        }
        throw new IllegalArgumentException("Couldn't do regression");
    }

    /**
     * returns the parameters B to minimize
     *  sum_i=1^n w_i * (y_i - (X_1 \dot B) ^ 2
     * 
     * Most of the matrix operations are unrolled but I use Colt
     * for the matrix inversion.
     */
    public double[] fitAnyVar(int first, int last) {
        /* ti = to invert = X^T * w * X */
        double ti[][] = new double[COLUMNS][COLUMNS];
        for (int i = 0; i < COLUMNS; i++) {
            XTWY[i] = 0;
            for (int j = 0; j < COLUMNS; j++) {
                ti[i][j] = 0;
            }
        }
        for (int i = first; i <= last; i++) {
            double wy = WY[i];
            for (int j = 0; j < COLUMNS; j++) {
                XTWY[j] += wy * X[i][j];
            }
        }
        for (int k = first; k <= last; k++) {
            for (int i = 0; i < COLUMNS; i++) {
                for (int j = 0; j <= i; j++) {
                    ti[i][j] += weights[k] * X[k][i] * X[k][j];
                }
            }            
        }
        for (int i = 0; i <= COLUMNS; i++) {
            for (int j = 0; j < COLUMNS; j++) {
                if (j > i) {
                    ti[i][j] = ti[j][i];
                }
            }
        }
        /* try the inversion in a loop.  If it fails because the matrix is 
           singular, add a little bit of noise and try again
        */
        double noisefactor = 0;
        Algebra alg = new Algebra();
        while (noisefactor < 10 * noiseIncrement) {
            try {
                DoubleMatrix2D toInvert = DoubleFactory2D.dense.make(ti);
                DoubleMatrix2D inverted = alg.inverse(toInvert);
                double output[] = new double[COLUMNS];
                for (int i = 0; i < COLUMNS; i++) {
                    output[i] = 0;
                    for (int j = 0; j < COLUMNS; j++) {
                        output[i] += XTWY[j] * inverted.getQuick(i,j);
                    }
                }
                return output;
            } catch (Exception e) {
                noisefactor += noiseIncrement;
                for (int i = 0; i < COLUMNS; i++) {
                    for (int j = 0; j < COLUMNS; j++) {
                        ti[i][j] += noisefactor * normal.nextDouble();
                    }
                }
            }
        }
        throw new IllegalArgumentException("Couldn't do regression");
    }

}
