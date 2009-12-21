/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

public class PaintableScaleChangedEvent {

	private PaintableScale scale;
	
	public PaintableScaleChangedEvent(PaintableScale sc) { 
		scale = sc;
	}
	
	public PaintableScale getScale() { return scale; }
}
