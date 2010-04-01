package edu.mit.csail.cgs.datasets.chipchip;

import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.utils.io.parsing.*;
import edu.mit.csail.cgs.utils.NotFoundException;

/* represents results from an agilent location array
   analysis along one chromosome. */

public interface ChipChipMLE extends GenericExperiment {
    public double getMax();
    /* start and stop are chromosomal positions */
    public double getMax(String chrom, int start, int stop) throws NotFoundException;
    /* start and stop are chromosomal positions.
       This sets the current window in the data to which
       all indexes i are relative.  Before the first call
       to window(), the window is the entire chromosome. */
    public void window(String chrom, int start, int stop, double minsize, double maxconf) throws NotFoundException;
    public double getSize(int i);
    public double getBindLL(int i);
    public double getNullLL(int i);
    public double getLogRat(int i);
    public double getConf(int i);

}
