/*
 * Author: tdanford
 * Date: Jun 9, 2008
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public interface AssignmentPredicate<X extends Region> {

	public boolean isValidAssignment(Region item, X event);
	public Mapper<Region,Region> assignmentZoneMapper();
	
	public static class Filter<Y extends Region> 
		implements edu.mit.csail.cgs.ewok.verbs.Filter<Y,Y> {
		
		private AssignmentPredicate<Y> predicate;
		private Region baseRegion;
		
		public Filter(AssignmentPredicate<Y> ap, Region r) {
			predicate = ap;
			baseRegion = r;
		}
		
		public Y execute(Y val) { 
			return predicate.isValidAssignment(baseRegion, val) ? val : null;
		}
	}
}


