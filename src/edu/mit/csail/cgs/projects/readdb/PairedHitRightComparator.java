package edu.mit.csail.cgs.projects.readdb;

import java.util.Comparator;

public class PairedHitRightComparator implements Comparator<PairedHit> {

    public boolean equals(Object o) {
        return (o instanceof PairedHitRightComparator);
    }

    public int compare(PairedHit a, PairedHit b) {
        if (a.rightChrom == b.rightChrom) {
            return a.rightPos - b.rightPos;
        } else {
            return a.rightChrom - b.rightChrom;
        }
    }
}