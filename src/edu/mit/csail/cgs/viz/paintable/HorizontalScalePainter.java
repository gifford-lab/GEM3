/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;

public class HorizontalScalePainter extends AbstractPaintable implements PaintableScaleListener {
	
	private PaintableScale scale;
	private boolean drawZero; 
	
	public HorizontalScalePainter(PaintableScale sc) { 
		scale = sc;
		scale.addPaintableScaleListener(this);
		drawZero = false;
	}
	
	public HorizontalScalePainter drawZero(boolean z) { 
		drawZero = z; 
		return this; 
	}

	public void paintItem(Graphics gg, int x1, int y1, int x2, int y2) {
		Graphics2D g = (Graphics2D)gg;
		g.setColor(Color.black);
		g.drawLine(x1, y2, x2, y2);
		
		int w = x2-x1, h = y2-y1;
		
		double min = scale.getMin(), max = scale.getMax();
		double range = max-min;
		
		double hash = Math.pow(10.0, (int)Math.floor(Math.log(range))-2.0);
		int ticks = (int)Math.floor(range/hash);
		
		//System.out.println(String.format("Hash %f -> %d ticks", hash, ticks));
		
		for(int i = 0; ticks >= 10; i++) { 
			hash *= i % 2 == 0 ? 5.0 : 2.0;
			ticks = (int)Math.floor(range/hash);
			//System.out.println(String.format("Hash %f -> %d ticks", hash, ticks));
		}
		
		double minTicks = Math.abs(min) / hash;
		double extra = minTicks - Math.floor(minTicks);
		double start = min + (hash-extra);
		//System.out.println(String.format("Min: %f, Extra: %f, Start: %f", min, extra, start));
		
		g.setColor(Color.black);
		if(drawZero && min < 0.0 && max > 0.0) { 
			double f = scale.fractionalOffset(0.0);
			int y = y2 - (int)Math.round(f * h);
			g.drawLine(x1, y, x2, y);
		}
		
		for(double hs = min + (hash-extra); hs <= max; hs += hash) { 
			double f = scale.fractionalOffset(hs);
			int x = x1 + (int)Math.round(f * w);
			
			g.setColor(Color.black);
			g.drawLine(x, y2, x, y2-5);
			
			g.translate(x, y2-5);
			g.rotate(-Math.PI/2.0);
			g.drawString(String.format("%.2f", hs), 0, 0);			
			g.rotate(Math.PI/2.0);
			g.translate(-x, -y2+5);
		}
	}

	public void paintableScaleChanged(PaintableScaleChangedEvent evt) {
		dispatchChangedEvent();
	}
}
