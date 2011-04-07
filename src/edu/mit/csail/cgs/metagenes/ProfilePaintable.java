/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
/**
 * 
 */
package edu.mit.csail.cgs.metagenes;

import java.awt.*;

import edu.mit.csail.cgs.viz.paintable.*;

/**
 * Wraps a Profile object, and paints it in a window region in the simplest possible way 
 * (connected lines between points).  
 * 
 * @author tdanford
 *
 */
public class ProfilePaintable extends AbstractPaintable implements ProfileListener {
	
	private PaintableScale scale;
	private BinningParameters params;
	private Profile profile;
	private Color col=Color.blue;
	private String style="Line";
	private boolean useNormalized=true;
	
	private int[] xs, ys;
	
	public ProfilePaintable(PaintableScale sc, Profile p) { 
		scale = sc;
		profile = p;
		profile.addProfileListener(this);
		params = profile.getBinningParameters();
		xs = new int[params.getNumBins()];
		ys = new int[params.getNumBins()];
	}
	
	public void profileChanged(ProfileEvent p) { 
		dispatchChangedEvent();
	}
	
	public Color getColor(){return col;}
	public void setColor(Color c){col=c;}
	public void setStyle(String s){style=s;}
	public void normalized(boolean n){useNormalized=n;}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke((float)2.5));
		
		int binPix = w / (params.getNumBins()+1);
		
		// Make sure that the profile doesn't change out from underneath us... 
		synchronized(profile) { 
			if(profile.min()<scale.getMin())
				scale.setScale(profile.min(), scale.getMax());
			if(profile.max()>scale.getMax())
				scale.setScale(scale.getMin(), profile.max());
			
			for(int i = 0; i < params.getNumBins(); i++) { 
				int x = x1 + (i+1)*binPix;
				double yf= scale.fractionalOffset(profile.value(i));
				int y = y2 - (int)Math.round(yf * (double)h);
				xs[i] = x; ys[i] = y;
			}
		
			//g2.setColor(Color.white);
			//g2.fillRect(x1, y1, w, h);
			g2.setColor(col);
			
			if(style.equals("Histo")){
				//Histogram style
				for(int i = 0; i < xs.length; i++) {
					g2.fillRect(xs[i],ys[i], binPix, h-ys[i]);
				}
			}else{
				//Line Graph style
				for(int i = 1; i < xs.length; i++) { 
					g2.drawLine(xs[i-1], ys[i-1], xs[i], ys[i]);
				}
				int rad = 2;
				int diam = rad*2;
				for(int i = 0; i < xs.length; i++) {
					g2.setColor(Color.white);
					g2.fillOval(xs[i]-rad, ys[i]-rad, diam, diam);
					g2.setColor(col);
					g2.drawOval(xs[i]-rad, ys[i]-rad, diam, diam);
				}
			}
		}
		g2.setStroke(oldStroke);
	}

}
