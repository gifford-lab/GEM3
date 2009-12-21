/*
 * Author: tdanford
 * Date: Sep 17, 2008
 */
package edu.mit.csail.cgs.viz.colors;

import java.awt.Color;

public abstract class Coloring {
	
	public static Color scale(Color min, Color max, double lambda) { 
		int r1 = min.getRed(), r2 = max.getRed();
		int g1 = min.getGreen(), g2 = max.getGreen();
		int b1 = min.getBlue(), b2 = max.getBlue();
		int a1 = min.getAlpha(), a2 = max.getAlpha();
		double l1 = 1.0 - lambda, l2 = lambda;
		int r = (int)Math.floor((double)r1 * l1 + (double)r2 * l2);
		int g = (int)Math.floor((double)g1 * l1 + (double)g2 * l2);
		int b = (int)Math.floor((double)b1 * l1 + (double)b2 * l2);
		int a = (int)Math.floor((double)a1 * l1 + (double)a2 * l2);
		return new Color(r, g, b, a);
	}
	
	public static Color opaque(Color c) { 
		int r = c.getRed(), g = c.getGreen(), b = c.getBlue();
		int a = 255;
		return new Color(r, g, b, a);
	}

	public static Color brighten(Color c) { 
		return c.brighter();
	}
	
	public static Color darken(Color c) { 
		return c.darker();
	}
	
	public static Color clearer(Color c) { 
		int[] rgba = rgba(c);
		rgba[3] /= 2;
		return asColor(rgba);
	}
	
	public static Color thicker(Color c) { 
		int[] rgba = rgba(c);
		rgba[3] = Math.min(255, rgba[3]*2);
		return asColor(rgba);
	}
	
	public static Color asColor(int[] rgb) { 
		if(rgb.length == 3) { 
			return new Color(rgb[0], rgb[1], rgb[2]);
		} else if (rgb.length == 4) { 
			return new Color(rgb[0], rgb[1], rgb[2], rgb[3]);			
		} else { 
			throw new IllegalArgumentException("array must be of length 3 or 4");
		}
	}
	
	public static int[] rgb(Color c) { 
		return new int[] { c.getRed(), c.getGreen(), c.getBlue() };
	}

	public static int[] rgba(Color c) { 
		return new int[] { c.getRed(), c.getGreen(), c.getBlue(), c.getAlpha() };
	}
}
