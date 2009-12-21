/*
 * Author: tdanford
 * Date: Sep 25, 2008
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

public class BitMasker {
	
	public static int[] bitMasks;
	
	static { 
		bitMasks = new int[32];
		int baseMask = 0x00000001;
		for(int i = 0; i < bitMasks.length; i++) {
			int mask = baseMask << i;
			bitMasks[i] = mask;
			//System.out.println(String.format("%d:  \t%08x", i, bitMasks[i]));
		}
	}
	
	public static int packBytes(byte[] b, int off) { 
		if(off < 0 || off+4 > b.length) throw new IllegalArgumentException();
		int blank = 0;
		blank += (((int)b[0]) << 24);
		blank += (((int)b[1]) << 16);
		blank += (((int)b[2]) << 8);
		blank += ((int)b[3]);
		return blank;
	}
	
	public static int setBit(int bits, int i, boolean value) { 
		if(value) { 
			return bits | bitMasks[i];
		} else { 
			return bits & (~bitMasks[i]);
		}
	}

	public static boolean hasAnyBit(int bits) { 
		return bits != 0;
	}
	
	public static boolean hasAllBits(int bits) { 
		return bits == 0xFFFFFFFF;
	}
	
	public static boolean hasBit(int bits, int i) {
		//System.out.println(String.format("%08x (@ %d : %08x) -> %08x", 
		//bits, i, bitMasks[i], (bits&bitMasks[i])));
		return (bits & bitMasks[i]) != 0;
	}
	
	public static Set<Integer> findBits(int bits) { 
		TreeSet<Integer> bitset = new TreeSet<Integer>();
		for(int i = 0; i < bitMasks.length; i++) { 
			if(hasBit(bits, i)) { 
				bitset.add(i);
			}
		}
		return bitset;
	}
	
	public static int countBits(int bits, boolean value) {
		int count = 0;
		for(int i = 0; i < bitMasks.length; i++) { 
			if(hasBit(bits, i) == value) { 
				count += 1;
			}
		}
		return count;
	}
	
	public static int findFirst(int bits, boolean value) { 
		for(int i = 0; i < bitMasks.length; i++) { 
			if(hasBit(bits, i) == value) { return i; }
		}
		return -1;
	}

	public static int findLast(int bits, boolean value) { 
		for(int i = bitMasks.length-1; i>= 0; i--) {  
			if(hasBit(bits, i) == value) { return i; }
		}
		return -1;
	}
}
