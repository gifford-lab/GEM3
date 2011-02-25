package edu.mit.csail.cgs.utils;

import java.util.List;
import java.util.ArrayList;

/**
 * Performs calculations for the UCSC bin indexing scheme.
 * See http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
 */

public class UCSCBins {

    public static final int[] binOffsets = {512+64+8+1, 64+8+1, 8+1, 1, 0};
    public static final int[] binSizes = {512*1024*1024, 64*1024*1024, 8*1024*1024, 1024*1024, 128*1024};
    public static final int binFirstShift = 17;       /* How much to shift to get to finest bin. */
    public static final int binNextShift = 3;         /* How much to shift to get to next larger bin. */

    /* Given start,end in chromosome coordinates assign it
     * a bin.   There's a bin for each 128k segment, for each
     * 1M segment, for each 8M segment, for each 64M segment,
     * and for each chromosome (which is assumed to be less than
     * 512M.)  A range goes into the smallest bin it will fit in. */
    public static int rangeToBin(int start, int end) throws IllegalArgumentException {        
        int startBin = start, endBin = end-1, i;
        if (start < 0) { throw new IllegalArgumentException("Start < 0");}
        if (end < 0) { throw new IllegalArgumentException("End < 0");}
        if (start > binSizes[0]) { throw new IllegalArgumentException("Start > 512M");}
        if (end > binSizes[0]) { throw new IllegalArgumentException("End > 512M");}


        startBin >>= binFirstShift;
        endBin >>= binFirstShift;
        for (i=0; i < binOffsets.length; ++i) {
            if (startBin == endBin) {
                return binOffsets[i] + startBin;
            }
            startBin >>= binNextShift;
            endBin >>= binNextShift;
        }
        throw new IllegalArgumentException(String.format("start %d, end %d out of range in findBin (max is 512M)", start, end));
    }
    /* returns the set of bins that cover a specified range */
    public static List<Integer> rangeToBins(int start, int end) throws IllegalArgumentException {
        if (start < 0) { throw new IllegalArgumentException("Start < 0");}
        if (end < 0) { throw new IllegalArgumentException("End < 0");}
        if (start > binSizes[0]) { throw new IllegalArgumentException("Start > 512M");}
        if (end > binSizes[0]) { throw new IllegalArgumentException("End > 512M");}
        List<Integer> output = new ArrayList<Integer>();
        output.add(0);
        for (int level = 1; level < binSizes.length; level++) {
            int binSize = binSizes[level];
            int first = start / binSize;
            int last = end / binSize;
            
            for (int i = first; i <= last; i++) {
                output.add(i + binOffsets[binOffsets.length - level - 1]);
            }
        }           
        return output;
    }
    public static String commaJoin(List<Integer> ints) {
        if (ints.size() == 0) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        sb.append(ints.get(0));
        for (int i = 1; i < ints.size(); i++) {
            sb.append("," + ints.get(i));
        }
        return sb.toString();
    }
    
}