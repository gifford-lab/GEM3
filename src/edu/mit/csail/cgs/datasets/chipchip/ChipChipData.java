package edu.mit.csail.cgs.datasets.chipchip;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.io.parsing.*;
import edu.mit.csail.cgs.utils.*;

/* represents results from an agilent location array
   along one chromosome. */

public interface ChipChipData extends GenericExperiment {

    public double getMax();
    public double getMin();
    /* start and stop are chromosomal positions */
    public double getMax(String chrom, int start, int stop) throws NotFoundException;
    public double getMin(String chrom, int start, int stop) throws NotFoundException;
    /* start and stop are chromosomal positions.
       This sets the current window in the data to which
       all indexes i are relative.  Before the first call
       to window(), the window is the entire chromosome. */
    public void window(String chrom, int start, int stop, double minratio) throws NotFoundException;
    /* get IP value for index i, replicate j */
    public double getRatio(int i, int j);
    public double getIP(int i, int j);
    public double getWCE(int i, int j);
    public double getVar(int i, int j);
    public char getStrand(int i, int j);
    /* returns an integer identifying the experiment from which the specified datapoint came */
    public int getExptID(int i, int j);
}
