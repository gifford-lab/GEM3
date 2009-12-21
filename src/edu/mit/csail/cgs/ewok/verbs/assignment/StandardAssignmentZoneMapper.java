/*
 * Author: tdanford
 * Date: Jun 9, 2008
 */
/**
 * 
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

/**
 * @author tdanford
 *
 * Turns a Region into the zone (another Region) into which an event must fall 
 * if it is to be annotated *to* the original Region.  
 */
public class StandardAssignmentZoneMapper implements Mapper<Region,Region> {
	
	private int upstream, downstream, nonstranded;
	
	public StandardAssignmentZoneMapper(int up, int down, int non) { 
		upstream = up;
		downstream = down;
		nonstranded = non;
	}
	
	public int getUpstream() { return upstream; }
	public int getDownstream() { return downstream; }
	public int getNonStranded() { return nonstranded; }

	public Region execute(Region a) {
		int start = a.getStart(), end = a.getEnd();
		if(a instanceof StrandedRegion) { 
			StrandedRegion sa = (StrandedRegion)a;
			start = sa.getStrand()=='+' ? sa.getStart() - upstream : sa.getEnd() - downstream;
			end = sa.getStrand() == '+' ? sa.getStart() + downstream : sa.getEnd() + upstream;
		} else { 
			start = a.getStart()-nonstranded; 
			end = a.getEnd()+nonstranded;
		}
		return new Region(a.getGenome(), a.getChrom(), start, end);
	}

}
