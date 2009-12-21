/*
 * Author: tdanford
 * Date: Nov 21, 2008
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

public class PackedBitVector implements BitVector {
	
	private int[] bits;
	private int length, extra;
	
	public PackedBitVector(int length) { 
		bits = new int[(length/32)+1];
		this.length = length;
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = 0;
		}
		
		extra = bits.length*32-length;
	}
	
	public PackedBitVector(PackedBitVector v) { 
		bits = new int[v.bits.length];
		length = v.length;
		for(int i =0; i < bits.length; i++) { 
			bits[i] = v.bits[i];
		}
	}
	
	public int[] asIndices(boolean value) {
		int[] array = new int[countOnBits()];

		for(int i = 0, j = 0; i < length; i++) {
			if(isOn(i) == value) { 
				array[j++] = i;
			}
		}
		return array;
	}
	
	public void copy(PackedBitVector v) { 
		if(v.length() != length()) { 
			throw new IllegalArgumentException(String.valueOf(v.length()));
		}
		
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = v.bits[i];
		}		
	}
	
	public void and(PackedBitVector v) { 
		if(v.length() != length()) { 
			throw new IllegalArgumentException(String.valueOf(v.length()));
		}
		
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = bits[i] & v.bits[i];
		}
	}

	public void or(PackedBitVector v) { 
		if(v.length() != length()) { 
			throw new IllegalArgumentException(String.valueOf(v.length()));
		}
		
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = bits[i] | v.bits[i];
		}
	}

	public void not() { 
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = ~bits[i];
		}
		
		clearExtraBits();
	}
	
	private void clearExtraBits() { 
		// Make sure we don't flip the "hanging" bits at the end.
		int last = bits.length-1;
		for(int k = 1; k <= extra; k++) { 
			bits[last] = BitMasker.setBit(bits[last], 32-k, false);
		}		
	}

	public int countOnBits() {
		int count = 0;
		for(int i = 0; i < bits.length; i++) { 
			count += BitMasker.countBits(bits[i], true);
		}
		return count;
	}

	public boolean isOff(int index) {
		int idx = index / 32;
		int offset = index % 32;
		return !BitMasker.hasBit(bits[idx], offset);
	}

	public boolean isOn(int index) {
		int idx = index / 32;
		int offset = index % 32;
		return BitMasker.hasBit(bits[idx], offset);
	}

	public int length() {
		return length;
	}

	public String toHexString() {
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < bits.length; i++) { 
			sb.append(String.format("%08x", bits[i]));
		}
		return sb.toString();
	}
	
	public void turnAllOn() { 
		for(int i = 0; i < bits.length; i++) { 
			bits[i] = 0xffffffff;
		}
		
		clearExtraBits();
	}
	
	public void turnAllOff() { 
		for(int i = 0;i < bits.length; i++) { 
			bits[i] = 0;
		}
	}
	
	public void turnOnRange(int i1, int i2) { 
		for(int i = i1; i < i2; i++) { 
			turnOnBit(i);
		}
	}

	public void turnOffRange(int i1, int i2) { 
		for(int i = i1; i < i2; i++) { 
			turnOffBit(i);
		}
	}

	public void turnOffBit(int index) {
		int idx = index / 32;
		int offset = index % 32;
		bits[idx] = BitMasker.setBit(bits[idx], offset, false);
	}

	public void turnOnBit(int index) {
		int idx = index / 32;
		int offset = index % 32;
		bits[idx] = BitMasker.setBit(bits[idx], offset, true);
	}

	public int hashCode() { 
		int code = 17;
		for(int i = 0; i < bits.length; i++) { 
			code += bits[i];
		}
		code *= 37;
		return code;
	}
	
	public boolean equal(Object o) { 
		if(!(o instanceof PackedBitVector)) { return false; }
		PackedBitVector v = (PackedBitVector)o;
		if(bits.length != v.bits.length) { return false; }
		if(length != v.length) { return false; }
		for(int i = 0; i < bits.length; i++) { 
			if(bits[i] != v.bits[i]) { 
				return false;
			}
		}
		return true;
	}
	
	public String toString() { 
		return String.format("%d: %s", length, toHexString());
	}
}
