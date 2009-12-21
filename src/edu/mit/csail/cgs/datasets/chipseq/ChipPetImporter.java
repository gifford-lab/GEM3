/*
 * Created on Sep 21, 2007
 */
package edu.mit.csail.cgs.datasets.chipseq;

import java.util.*;
import edu.mit.csail.cgs.datasets.chippet.*;

/**
 * @author tdanford
 *
 * Converts entries from the old chippet schema into the chipseq schema of the future -- 
 * this simultaneously creates 
 * (1) an experiment
 * (2) a set of reads for that experiment
 * (3) an alignment for those reads to a particular genome, for that experiment
 * (4) entries in the chipseqhits table.
 */
public class ChipPetImporter {

    private ChipPetLoader cpLoader;
    private ChipSeqLoader csLoader;
    
    public ChipPetImporter() { 
        
    }
}
