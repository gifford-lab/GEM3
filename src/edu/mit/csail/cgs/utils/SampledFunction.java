package edu.mit.csail.cgs.utils;

import java.lang.*;
import java.io.*;
import java.util.*;

public class SampledFunction { 
    
    private double[] fValues;
    private double[] fPoints;
    private int[] fPixelValues;
    private int fPosSpace, fNegSpace;
    private double fMax, fMin, fScale;
    private boolean fRecalcScale;
    
    private RealValuedFunction fFunc;

    public SampledFunction(RealValuedFunction f, double[] pts, 
			   int posSpace, int negSpace,
			   double scale) { 
	fFunc = f;
	fRecalcScale = false;
	fScale = scale;
	resampleFunction(pts, posSpace, negSpace);
    }

    public String getName() { return fFunc.getName(); }
    public double getScale() { return fScale; }

    public void setScale(double v) { 
	fScale = v; 
	fRecalcScale = false;
    }

    public void setRecalcScale(boolean v) { fRecalcScale = v; }

    public void resampleFunction(double[] pts, int posSpace, int negSpace) { 

	//System.out.println("Resampling Function: " + 
	//		   posSpace + "," + negSpace);

	fPosSpace = posSpace;
	fNegSpace = negSpace;
	fPoints = (double[])pts.clone();
	fValues = new double[fPoints.length];
	fPixelValues = new int[fPoints.length];
	fMax = -Double.MAX_VALUE;
	fMin = Double.MAX_VALUE;

	//System.out.print("\t["); System.out.flush();
	for(int i = 0; i < fValues.length; i++) { 
	    try { 

		fValues[i] = fFunc.eval(fPoints[i]);
		if(fValues[i] > fMax) { fMax = fValues[i]; }
		if(fValues[i] < fMin) { fMin = fValues[i]; }

	    } catch(IllegalArgumentException iae) { 
		fValues[i] = 0.0;  // but doesn't affect min/max
	    }
	    //System.out.print("."); System.out.flush();
	}
	//System.out.println("] (" + fMax + "," + fMin + ")");

	double posScale = Double.MAX_VALUE;
	double negScale = Double.MAX_VALUE;
	if(fMax > 0.0) { 
	    posScale = (double)fPosSpace / fMax;
	} 
	
	if(fMin < 0.0) { 
	    negScale = (double)(fNegSpace) / (-fMin);
	}
	
	double naturalScale = Double.MAX_VALUE;
	if(posScale < naturalScale) { naturalScale = posScale; }
	if(negScale < naturalScale) { naturalScale = negScale; }

	//System.out.println("\tScale (" + posScale + " / " + negScale + 
	//		   ") --> " + fScale);

	//System.out.print("\tPixels: ["); System.out.flush();
	//System.out.println("Resampled Function Scale: " + fScale);
	for(int i = 0; i < fValues.length; i++) { 
	    fPixelValues[i] = 
		(int)Math.round(fValues[i] * naturalScale * fScale);
	    //System.out.print("."); System.out.flush();
	}
	//System.out.println("]");
    }

    public int size() { return fPoints.length; }
    public double getValue(int idx) { return fValues[idx]; }
    public int getPixelValue(int idx) { return fPixelValues[idx]; }
    public double getPoint(int idx) { return fPoints[idx]; }
}

