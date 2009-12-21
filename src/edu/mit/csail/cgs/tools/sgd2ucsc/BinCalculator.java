/*
 * Created on Jan 26, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.sgd2ucsc;

public class BinCalculator {

    private int binFirstShift; 
    private int binNextShift;
    private int[] binOffsets;
    
    public BinCalculator() { 
        binFirstShift = 17;
        binNextShift = 3;
        binOffsets = new int[5];
        binOffsets[0] = 512+64+8+1;
        binOffsets[1] = 64+8+1;
        binOffsets[2] = 8+1;
        binOffsets[3] = 1;
        binOffsets[4] = 0;
    }
    
    public BinCalculator(int fs, int ns, int[] bo) {
        binFirstShift = fs;
        binNextShift = ns;
        binOffsets = (int[])bo.clone();
    }
 
    public int getBinFromRange(int start, int end) { 
        int startbin = start;
        int endbin = end-1;
        startbin = startbin >> binFirstShift;
        endbin = endbin >> binFirstShift;
        
        for(int i = 0; i < binOffsets.length; i++) { 
            if(startbin == endbin) { 
                return binOffsets[i] + startbin;
            }
            startbin = startbin >> binNextShift;
            endbin = endbin >> binNextShift;
        }
        
        System.err.println("Warning: Can't get bin for " + 
                start + "," + end + " --> " + startbin + "," + endbin);
        return 0;
    }
}
