package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.utils.*;

/** superclass for any dataset that presents values at various
 *  genomic coordinates 
 */

public interface GenericExperiment {
    
    public String getName();
    public String getVersion();
    
    /**
     * set the current data window.   
     */
    public void window(String chrom, int start, int stop) throws NotFoundException ;
    /** returns the number of data indices in this window */
    public int getCount();
    /** returns the number of replicates at index i*/
    public int getReplicates(int i);
    /** returns the value at index i, replicate j */
    public double getValue(int i, int j);
    /** returns the genomic coordinate of the observation 
        at index i*/
    public int getPos(int i);
    /** convert a base to the nearest index.  Returns -1
        if the base pair position isn't in this window */
    public int baseToIndex(int b);
    public String getWindowChrom();
    public int getWindowStart();
    public int getWindowStop();
    public double getMax(String chrom, int start, int stop) throws NotFoundException;
    public double getMin(String chrom, int start, int stop) throws NotFoundException;
}
