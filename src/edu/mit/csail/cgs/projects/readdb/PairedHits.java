package edu.mit.csail.cgs.projects.readdb;

import java.io.*;

/**
 * Represents a list of sorted reads on disk
 */
public class PairedHits extends Hits {

    /** 
     * true if the positions (and chrom) for these PairedHits are
     * the left side, meaning that the other chrom and position are for the right side.
     * If false, then the positions are for the right side and the other chrom and
     * position are for the left side.
     */
    private boolean isLeft;
    


}