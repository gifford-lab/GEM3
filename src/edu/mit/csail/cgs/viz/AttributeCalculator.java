/*
 * Created on Apr 3, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package edu.mit.csail.cgs.viz;

/**
 * @author Lindy Briggs
 *
 * AttributeCalculator is an interface, originally created by Lindy,
 * that allows us to get dynamic "attributes" (line widths, point widths, etc)
 * based on the *size* of the window that we're drawing into.  This allows for 
 * a more natural scaling process, and doesn't require us to write the same 
 * sizing code into each of our painters (or worse, have them all use a fixed 
 * style at all sizes).
 */
public interface AttributeCalculator {
	
	public int getPointWidth(int w, int h, double scale);
	public int getLineWidth(int w, int h, double scale);
	public int getFontSize(int w, int h, double scale);
	
	public void setPointWidth(int w);
	public void setLineWidth(int w);
	public void setFontSize(int size);

}
