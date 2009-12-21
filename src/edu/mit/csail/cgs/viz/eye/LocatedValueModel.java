/*
 * Author: tdanford
 * Date: Aug 7, 2009
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.Color;

import edu.mit.csail.cgs.utils.models.Model;

public class LocatedValueModel extends Model {
	
	public Integer location;
	public Double value; 
	public Color color;
	
	public LocatedValueModel() {}
	
	public LocatedValueModel(int loc, double val, Color c) { 
		location = loc; value = val;
		color = c;
	}
}