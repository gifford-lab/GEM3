package edu.mit.csail.cgs.viz.expression;

import java.awt.Color;

public class CombinationColorMap implements IntensityColorMap {
	
	private double min, max, split;
	private IntensityColorMap lower, upper;
	
	public CombinationColorMap(double min, double split, 
			double max, IntensityColorMap lower, IntensityColorMap upper) { 
		this.min = min; this.max = max;
		this.split = split;
		this.lower = lower;
		this.upper = upper;
		
		if(split < min || split > max) { 
			throw new IllegalArgumentException();
		}
	}

	public Color getColor(double v) {
		if(v >= min && v <= max) { 
			if(v <= split) { 
				double f = (v - min) / (split-min);
				double nv = lower.getMinValue() + 
					(f * (lower.getMaxValue()-lower.getMinValue()));
				return lower.getColor(nv);
			} else { 
				double f = (v - split) / (max-split);
				double nv = upper.getMinValue() + 
					(f * (upper.getMaxValue()-upper.getMinValue()));
				return upper.getColor(nv);				
			}
		} else { 
			return Color.white; 
		}
	}

	public double getMaxValue() {
		return max;
	}

	public double getMinValue() {
		return min;
	}

}
