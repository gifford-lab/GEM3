package edu.mit.csail.cgs.projects.readdb;

import java.util.Comparator;

public class PairedHitLeftComparator implements Comparator<PairedHit> {

    public boolean equals(Object o) {
        return (o instanceof PairedHitLeftComparator);
    }

    public int compare(PairedHit a, PairedHit b) {
        int result;
        if (a.leftChrom == b.leftChrom) {
            result = a.leftPos - b.leftPos;
        } else {
            result = a.leftChrom - b.leftChrom;
        }
        if (result == 0) {
            if (a.rightChrom == b.rightChrom) {
                result = a.rightPos - b.rightPos;
            } else {
                result = a.rightChrom - b.rightChrom;
            }
        }
        if (result == 0) {
            result = a.leftLength - b.leftLength;
        }
        if (result == 0) {
            result = a.rightLength - b.rightLength;
        }
        return result;
    }
}