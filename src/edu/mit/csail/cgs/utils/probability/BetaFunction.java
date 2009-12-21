package edu.mit.csail.cgs.utils.probability;

import java.lang.*;
import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.RealValuedFunction;

public class BetaFunction implements RealValuedFunction { 

    private static double cof[];
    
    static {
        cof = new double[6];
        cof[0] = 76.18009172947146;
        cof[1] = -86.5053203941677;
        cof[2] = 24.01409824083091;
        cof[3] = -1.231739572450155;
        cof[4] = 0.1208650973866179e-2;
        cof[5] = -0.5395239384953e-5;
    }   
    
    private double fAlpha, fBeta;
    private String fName;
    
    public BetaFunction(double a, double b) { 
        fAlpha = a; fBeta = b;
        if(fAlpha <= 0.0 || fBeta <= 0.0) { 
            throw new IllegalArgumentException();
        }
        fName = "Beta(" + a + "," + b + ")";
    }
    
    public double eval(double input) { 
        if(input < 0 || input > 1.0) { 
            throw new IllegalArgumentException(); 
        }
        
        double f1 = Math.pow(input, (fAlpha - 1.0));
        double f2 = Math.pow(1.0 - input, fBeta - 1.0);
        double c = Math.exp(gammln(fAlpha + fBeta) - 
                (gammln(fAlpha) + gammln(fBeta)));
        
        return c * f1 * f2;
    }
    
    private double gammln(double xx) { 
        double x, y, tmp, ser;
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * Math.log(tmp);
        ser = 1.000000000190015;
        for(j = 0; j <= 5; j++) { 
            ser += cof[j] / ++y;
        }
        return -tmp + Math.log(2.5066282746310005 * ser / x);
    }
    
    public String getName() { 
        return fName;
    }
}

