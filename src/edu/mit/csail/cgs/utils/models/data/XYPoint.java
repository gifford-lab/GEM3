/*
 * Author: tdanford
 * Date: Aug 28, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import edu.mit.csail.cgs.utils.models.Model;


public class XYPoint extends Model {
	public Double x, y;
	
	public XYPoint() {}
	public XYPoint(Double _x, Double _y) { x = _x; y = _y; }
}
