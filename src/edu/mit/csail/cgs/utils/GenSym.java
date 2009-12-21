/*
 * Author: tdanford
 * Date: May 27, 2008
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * Generates unique string symbols. 
 * 
 * @author tdanford
 */
public class GenSym {

	private String prefix;
	private long index;
	
	public GenSym(String pre) { 
		prefix = pre;
		index = (long)0;
	}
	
	public synchronized String generateSymbol() { 
		String sym = String.format("%s%d", prefix, index);
		index += 1;
		return sym;
	}
}
