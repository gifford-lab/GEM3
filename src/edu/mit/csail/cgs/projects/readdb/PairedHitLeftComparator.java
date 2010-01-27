package edu.mit.csail.cgs.projects.readdb;

import java.util.Comparator;

public class PairedHitLeftComparator implements Comparator<PairedHit> {

    public boolean equals(Object o) {
        return (o instanceof PairedHitLeftComparator);
    }

    public int compare(PairedHit a, PairedHit b) {
        if (a.leftChrom == b.leftChrom) {
            return a.leftPos - b.leftPos;
        } else {
            return a.leftChrom - b.leftChrom;
        }
    }
}