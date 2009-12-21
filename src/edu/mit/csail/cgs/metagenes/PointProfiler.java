/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.mit.csail.cgs.metagenes;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.ewok.verbs.Filter;

public interface PointProfiler<PointClass extends Point, ProfileClass extends Profile> extends Filter<PointClass,ProfileClass> { 
	public BinningParameters getBinningParameters();
	public void cleanup();
}
