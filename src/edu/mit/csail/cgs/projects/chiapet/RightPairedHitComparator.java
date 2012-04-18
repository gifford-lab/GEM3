package edu.mit.csail.cgs.projects.chiapet;

import java.util.Comparator;

import edu.mit.csail.cgs.projects.readdb.PairedHit;

public class RightPairedHitComparator implements Comparator<PairedHit> {
	
	public int compare(PairedHit arg0, PairedHit arg1) {
		int tor;
		if (arg0.rightChrom == arg1.rightChrom) {
			tor = arg0.rightPos - arg1.rightPos;
		} else {
			tor = arg0.rightChrom - arg1.rightChrom;
		}
		if (tor==0) {
			if (arg0.leftChrom==arg1.leftChrom) {
				return arg0.leftPos-arg1.leftPos;
			} else {
				return arg0.leftChrom - arg1.leftChrom;
			}
		} else {
			return tor;
		}
	}

}
