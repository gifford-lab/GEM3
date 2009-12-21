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
 * @author tdanford
 *
 */
public class ProfileLinePaintable extends AbstractPaintable implements ProfileListener {
	
	private PaintableScale scale;
	private BinningParameters params;
	private Profile profile;
	private Color col=Color.blue;
	private boolean quantized=false;
	private double[] quanta;
	
	public ProfileLinePaintable(PaintableScale sc, Profile p) { 
		scale = sc;
		profile = p;
		profile.addProfileListener(this);
		params = profile.getBinningParameters();
	}
	
	public Profile getProfile() { return profile; }
	public void setColor(Color c){col=c;}
	
	public void setQuanta(double [] q){
		if(q!=null){
			quantized=true;
			quanta=q;
		}
	}
	
	public void profileChanged(ProfileEvent p) { 
		dispatchChangedEvent();
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		Graphics2D g2 = (Graphics2D)g;
		
		int binPix = w / (params.getNumBins());

		for(int i = 0; i < params.getNumBins(); i++) { 
			int x = x1 + i*binPix;
			double value = profile.value(i);
			double yf;
			
			if(quantized){
				yf=0;
				for(int a=0; a<quanta.length; a++)
					if(value>quanta[a])
						yf=(double)a/((double)quanta.length-1);
			}else{
				yf = scale.fractionalOffset(value);
			}
			
			Color c = calcFracColor(col, yf);
			g2.setColor(c);
			g2.fillRect(x, y1, binPix, h);
		}
	}
	
	private Color calcFracColor(Color col, double v){
		Color c;
		Color maxColor = col;
		Color minColor = Color.white;
		
		double sVal = v>1 ? 1 : (v<0 ? 0 : v);
		int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
	    int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
	    int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
	    c = new Color(red, green, blue);
		return(c);
	}
	
}
