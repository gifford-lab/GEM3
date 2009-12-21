/*
 * Author: tdanford
 * Date: Jun 9, 2008
 */
/**
 * 
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

/**
 * @author tdanford
 *
 */
public class StandardAssignmentPredicate<X extends Region> implements AssignmentPredicate<X> {
	
	private int upstream, downstream, nonstranded;
	
	public StandardAssignmentPredicate(int up, int down, int ns) { 
		upstream = up;
		downstream = down;
		nonstranded = ns;
	}

	public Mapper assignmentZoneMapper() {
		return new StandardAssignmentZoneMapper(upstream, downstream, nonstranded);
	}
	
	private boolean contains(Region r, int pt) { 
		return r.getStart() <= pt && r.getEnd() >= pt;
	}

	public boolean isValidAssignment(Region item, X event) {
		if(!item.getChrom().equals(event.getChrom())) { return false; }
		if(item instanceof StrandedRegion) { 
			StrandedRegion sitem = (StrandedRegion)item;
			boolean strand = sitem.getStrand() == '+';
			if(strand) { 
				if(contains(event, sitem.getStart())) { return true; }
				if(event.getEnd() < item.getStart()) { 
					int dist = item.getStart()-event.getEnd();
					return dist <= upstream;
				} else { 
					int dist = event.getStart()-item.getStart();
					return dist <= downstream;
				}
			} else { 
				if(contains(event, sitem.getEnd())) { return true; }
				if(event.getEnd() < item.getEnd()) { 
					int dist = item.getEnd()-event.getEnd();
					return dist <= downstream;
				} else { 
					int dist = event.getStart()-item.getEnd();
					return dist <= upstream;
				}				
			}
		} else { 
			int dist = item.distance(event);
			return dist <= nonstranded;
		}
	}
}
