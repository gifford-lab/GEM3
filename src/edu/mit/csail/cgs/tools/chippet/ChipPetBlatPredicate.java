/*
 * Created on May 14, 2007
 */
package edu.mit.csail.cgs.tools.chippet;

import edu.mit.csail.cgs.utils.parsing.alignment.BlatPSLEntry;
import edu.mit.csail.cgs.utils.parsing.alignment.BlatPSLEntryPredicate;

public class ChipPetBlatPredicate implements BlatPSLEntryPredicate { 

    private int minMatches;
    private int minGap;

    public ChipPetBlatPredicate(int mm, int mg) { 
        minMatches = mm;
        minGap = mg;
    }

    public ChipPetBlatPredicate() { 
        minMatches = 30;
        minGap = 10;
    }

    public boolean acceptsEntry(BlatPSLEntry e) { 
        if(e.getMatch() < minMatches) { return false; }
        if(e.getQgapCount() > 0) { return false; }
        if(e.getTgapCount() != 1) { return false; }
        if(e.getTGapSize(0) < minGap) { return false; }
        return true;
    }
}
