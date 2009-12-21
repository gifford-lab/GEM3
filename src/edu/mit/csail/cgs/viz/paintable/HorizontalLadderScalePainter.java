/*
 * Author: tdanford
 * Date: Jul 18, 2008
 */
package edu.mit.csail.cgs.viz.paintable;

import java.awt.*;

public class HorizontalLadderScalePainter extends AbstractPaintable implements PaintableScaleListener {
	
	private PaintableScale scale;
	
	public HorizontalLadderScalePainter(PaintableScale sc) { 
		scale = sc;
		scale.addPaintableScaleListener(this);
	}
	
	public HorizontalLadderScalePainter drawZero(boolean z) { 
		return this; 
	}

	public void paintItem(Graphics gg, int x1, int y1, int x2, int y2) {
		Graphics2D g = (Graphics2D)gg;
		
		int w = x2-x1, h = y2-y1;
		int my = y1 + h/2;

		g.setColor(Color.black);
		g.drawLine(x1, my, x2, my);
		
		int min = (int)Math.floor(scale.getMin()), max = (int)Math.ceil(scale.getMax());
		int range = max-min;
		
		int hashSize = 1;
		int ticks = range/hashSize;
		
		//System.out.println(String.format("Hash %f -> %d ticks", hash, ticks));
		
		for(int i = 0; ticks >= 20; i++) { 
			hashSize *= i % 2 == 0 ? 5 : 2;
			ticks = range / hashSize;
		}
		
		int diff = min % hashSize;
		int extra = diff > 0 ? (hashSize - diff) : 0;
		
		g.setColor(Color.black);
		
		FontMetrics fm = g.getFontMetrics();
		Font oldFont = g.getFont();
		int lastx = x1 - 100;

		for(int hs = min + extra; hs <= max; hs += hashSize) { 
			double f = scale.fractionalOffset((double)hs);
			int x = x1 + (int)Math.round(f * w);
			
			g.setColor(Color.black);
			g.drawLine(x, my-5, x, my);
			
			String str = String.format("%d", hs);
			int strwidth = fm.charsWidth(str.toCharArray(), 0, str.length());
			int strheight = fm.getAscent();

			if(x - strwidth/2 > lastx) { 
				g.translate(x-strwidth/2, my + strheight + 2);
				g.drawString(str, 0, 0);			
				g.translate(-x+strwidth/2, -my - strheight - 2);
				lastx = x + strwidth/2;
			}
		}
		
		g.setFont(oldFont);
	}

	public void paintableScaleChanged(PaintableScaleChangedEvent evt) {
		dispatchChangedEvent();
	}
}
