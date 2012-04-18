package edu.mit.csail.cgs.projects.chiapet;

import java.util.Comparator;

import edu.mit.csail.cgs.projects.readdb.PairedHit;

public class LeftPairedHitComparator implements Comparator<PairedHit> {

	public int compare(PairedHit arg0, PairedHit arg1) {
		int tor;
		if (arg0.leftChrom==arg1.leftChrom) {
			tor = arg0.leftPos-arg1.leftPos;
		} else {
			tor = arg0.leftChrom - arg1.leftChrom;
		}
		if (tor==0) {
			if (arg0.rightChrom == arg1.rightChrom) {
				return arg0.rightPos - arg1.rightPos;
			} else {
				return arg0.rightChrom - arg1.rightChrom;
			}
		} else {
			return tor;
		}
	}

}
