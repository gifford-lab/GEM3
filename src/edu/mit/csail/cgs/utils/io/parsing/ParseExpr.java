/*
 * Created on May 23, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;

/**
 * @author tdanford
 */
public interface ParseExpr { 
    public String getExptName(int i);
    public Map<String,Double> getWholeExpt(int i);
    public boolean exptHasORF(int i, String orf);
    public double getValue(int i, String orf);
	public Double[] getValues(int start, int end, String orf);
    public int numExpts();
    public int numORFs(int i);
    public double getMax(int i);
    public double getMin(int i);
}
