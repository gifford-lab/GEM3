/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.util.*;
import java.awt.*;

public class PaintableScale {
	
	private double min, max;
	private LinkedList<PaintableScaleListener> listeners; 
	
	public PaintableScale(double m1, double m2) { 
		min = m1; 
		max = m2;
		listeners = new LinkedList<PaintableScaleListener>();
	}
	
	public String toString() { 
		return String.format("Scale: %.3f-%.3f", min, max);
	}
	
	public void addPaintableScaleListener(PaintableScaleListener psl) { 
		listeners.addLast(psl);
	}
	
	public void removePaintableScaleListener(PaintableScaleListener psl) { 
		listeners.remove(psl);
	}

	public double getMin() { return min; }
	public double getMax() { return max; }
	
	public double getRange() { return max - min; }
	public double fractionalOffset(double value) { return (value-min)/getRange(); }
	
	public void setScale(double m1, double m2) { 
		if(m1 > m2) { throw new IllegalArgumentException(); }
		min = m1; 
		max = m2; 
		dispatchPaintableChangedEvent();
	}
	
	public void updateScale(double newPoint) { 
		if(newPoint < min){ 
			min = newPoint;
			dispatchPaintableChangedEvent();
		}else if(newPoint>max){
			max = Math.max(newPoint, max);
			dispatchPaintableChangedEvent();
		}
	}
	
	private void dispatchPaintableChangedEvent() { 
		PaintableScaleChangedEvent evt = new PaintableScaleChangedEvent(this);
		for(PaintableScaleListener list : listeners) { 
			list.paintableScaleChanged(evt);
		}
	}
	
	private static double pow(double base, int exp) { 
		double x = 1.0;
		for(int i = 0; i < exp; i++) { x *= base; }
		for(int i = 0; i > exp; i--) { x /= base; }
		return x;
	}

	public double[] findHashMarks(int maxTicks) { 
		double range = max-min;
		double exp = Math.log(range) / Math.log(10.0);
		double hash = pow(10.0, (int)Math.floor(exp));
		int ticks = (int)Math.floor(range/hash);
		
		double[] mults = new double[] { 5.0, 2.0 };
		
		for(int i = 0; ticks > maxTicks; i += 1) { 
			hash *= mults[i % mults.length];
			ticks = (int)Math.floor(range / hash);
		}
		
		double lower = Math.floor(min/hash) * hash;
		if(lower < min) { lower += hash; }
		
		double upper = Math.ceil(max/hash) * hash;
		if(upper > max) { upper -= hash; }

		int count = (int)Math.round((upper-lower)/hash);
		double[] array = new double[count];
		
		double h = lower;
		for(int i = 0; i < array.length; i++, h += hash) { 
			array[i] = h;
		}
		
		return array;
	}
}
