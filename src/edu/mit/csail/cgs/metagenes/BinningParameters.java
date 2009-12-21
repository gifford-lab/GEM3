/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.mit.csail.cgs.metagenes;

import java.util.*;
import java.util.regex.*;

public class BinningParameters {
	
	private int windowSize, bins, binSize;	
	
	public BinningParameters(int window, int nbins) { 
		windowSize = window;
		bins = nbins;
		binSize = windowSize / bins;
		if(windowSize % binSize != 0) { 
			throw new IllegalArgumentException(String.format("Window size %d " +
					"must be a multiple of bin count %d", windowSize, bins));
		}
	}
	
	public BinningParameters(String ps) { 
		Pattern p = Pattern.compile("(\\d+)/(\\d+)");
		Matcher m = p.matcher(ps);
		if(m.find()) { 
			windowSize = Integer.parseInt(m.group(1));
			bins = Integer.parseInt(m.group(2));
			binSize = windowSize / bins;

			if(windowSize % binSize != 0) { 
				throw new IllegalArgumentException(String.format("Window size %d " +
						"must be a multiple of bin count %d", windowSize, bins));
			}
		} else {
			throw new IllegalArgumentException(String.format("%s is in incorrect format", ps));
		}
	}
	
	public int getWindowSize() { return windowSize; }
	public int getNumBins() { return bins; }
	public int getBinSize() { return binSize; }
	
	public int findBin(int offset) { 
		return Math.max(0, Math.min(bins-1, offset/binSize));
	}
}
