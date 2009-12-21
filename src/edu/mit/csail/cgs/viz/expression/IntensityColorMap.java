package edu.mit.csail.cgs.viz.expression;

import java.awt.Color;

public interface IntensityColorMap {
	public double getMinValue();
	public double getMaxValue();
	public Color getColor(double v);
	
	public static class LinearColorMap implements IntensityColorMap {
		
		private double minValue, maxValue;
		private int minColor, maxColor;
		private int index;
		
		public LinearColorMap(int rgb, int minCol, int maxCol, double minVal, double maxVal) { 
			index = rgb;
			minColor = minCol;
			maxColor = maxCol;
			minValue = minVal;
			maxValue = maxVal;
		}
		
		public double getMinValue() { return minValue; }
		public double getMaxValue() { return maxValue; }
		
		public int getMinColor() { return minColor; }
		public int getMaxColor() { return maxColor; }
		public int getIndex() { return index; }
		
		public Color getColor(double v) { 
			int[] carray = new int[3];
			carray[0] = carray[1] = carray[2] = 0;
			
			double f = (v - minValue) / (maxValue - minValue);
			if(f <= 0.0) { f = 0.0; }
			if(f > 1.0) { f = 1.0; }
			carray[index] = minColor + (int)Math.floor(f * (double)(maxColor-minColor));
			
			return new Color(carray[0], carray[1], carray[2]);
		}
	}
	
	public static class RedMap extends LinearColorMap { 
		public RedMap() { 
			super(0, 0, 255, 0.0, 1.0);
		}
		
		public RedMap(double min, double max) { 
			super(0, 0, 255, min, max);
		}
	}
	
	public static class GreenMap extends LinearColorMap { 
		public GreenMap() { 
			super(0, 0, 255, 0.0, 1.0);
		}

		public GreenMap(double min, double max) { 
			super(1, 0, 255, min, max);
		}
	}
	
	public static class BlueMap extends LinearColorMap { 
		public BlueMap() { 
			super(0, 0, 255, 0.0, 1.0);
		}

		public BlueMap(double min, double max) { 
			super(2, 0, 255, min, max);
		}
	}
}
