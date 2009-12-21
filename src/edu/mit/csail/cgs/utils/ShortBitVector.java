/*
 * Created on Jun 20, 2005
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

/**
 * @author tdanford
 */
public class ShortBitVector implements Comparable<ShortBitVector>, BitVector {
    
    public static void main(String[] args) { 
        String b = args[0];
        char[] array = b.toCharArray();
        char m = '1';
        short v = ON;
        ShortBitVector bv = new ShortBitVector(array, m, v);
        String expandString = args[1];
        System.out.println(expandString + " --> " + bv.expandString(expandString, "-"));
    }
    
    public static short ON = (short)1;
    public static short OFF = (short)0;
    
    private static char[] hexChars = {'0', '1', '2', '3', '4', '5', '6', '7', 
        '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'
    };
    
    private short[] array;

    /**
     * Creates a blank bit-vector
     * @param length The length of the blank vector.
     */
    public ShortBitVector(int length) {
        array = new short[length];
        for(int i = 0; i < length; i++) { 
            array[i] = OFF;
        }
    }
    
    /**
     * Creates a blank bit-vector
     * @param length The length of the blank vector.
     * @param value The value of each bit in the new vector.
     */
    public ShortBitVector(int length, short value) {
        if(value != ON && value != OFF) { throw new IllegalArgumentException(); }
        array = new short[length];
        for(int i = 0; i < length; i++) { 
            array[i] = value;
        }
    }
    
    /**
     * Creates a pre-set Bit-vector
     * @param _array The array of values for the vector: 0 values become OFF, 
     * everything else becomes ON.
     */
    public ShortBitVector(int[] _array) { 
        array = new short[_array.length];
        for(int i = 0; i < array.length; i++) { 
            if(_array[i] == 0) { 
                array[i] = OFF;
            } else { 
                array[i] = ON;
            }
        }
    }
    
    /**
     * Creates a BitVector from a raw string representation, with length equal to the
     * character array supplied.  If _array[i] == maker, then bit_vector[i] == value.
     * @param _array The character array, whose values dictate the contents of the bit vector.
     * @param marker The character that indicates which bits have 'value.'
     * @param value The value of all bits that correspond to "marker" character.
     */
    public ShortBitVector(char[] _array, char marker, short value) { 
        if(value != ON && value != OFF) { throw new IllegalArgumentException(); }
        array = new short[_array.length];
        short other_value = OFF;
        if(value == OFF) { other_value = ON; }
        for(int i = 0; i < array.length; i++) { 
            if(_array[i] == marker) { 
                array[i] = value;
            } else { 
                array[i] = other_value;
            }
        }
    }
    
    /**
     * Copy constructor.
     * @param bv The bit-vector to copy.
     */
    public ShortBitVector(ShortBitVector bv) { 
        array = new short[bv.array.length];
        for(int i = 0; i < array.length; i++) { array[i] = bv.array[i]; }
    }
    
    /**
     * Converts a hex string into a BitVector.
     * @param length The length of the bit-vector 
     * @param hArray The hex-character string.
     */
    public ShortBitVector(int length, char[] hArray) { 
        int starter = length % 4;
        array = new short[length];
        int hIndex = 0;
        int aIndex = 0;
        if(starter > 0) { 
            translateHex2Short(hArray[hIndex], array, aIndex, starter);
            hIndex += 1;
            aIndex += starter;
        }
        for(; aIndex < length; aIndex += 4, hIndex += 1) { 
            translateHex2Short(hArray[hIndex], array, aIndex, 4);
        }
    }
    
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#toHexString()
	 */
    public String toHexString() { 
        StringBuilder sb = new StringBuilder();
		int index = 0;
        int start = array.length % 4;
        if(start > 0) { 
            sb.append(translateShort2Hex(array, index, start));
			index += start;
        }
        for(; index < array.length; index += 4) { 
            sb.append(translateShort2Hex(array, index, 4));
        }
        return sb.toString();
    }
    
    /**
     * expandString(base, offInsert) takes a string, which should have a length equal to 
     * the number of ON bits in the vector.  It then returns the same string, but with 
     * the substring offInsert inserted repeatedly for each occurrence of an OFF bit.
     * For instance, the bit vector 1101101 will expand the string "aabbc" (with offInsert
     * set to "-") into the string "aa-bb-c".
     * 
     * @param base  The "base" string which will be inserted.
     * @param offInsert The substring which we will insert for each occurrence of an OFF
     * bit in the bit-vector.
     * @return Returns the expanded String.
     */
    public String expandString(String base, String offInsert) {
        if(countOnBits() != base.length()) { 
            throw new IllegalArgumentException("length " + base.length() + " != " + 
                    countOnBits());
        }
        StringBuilder sb = new StringBuilder(base);
        int bi = 0, index = 0;
        int offLength = offInsert.length();
        for(; index < array.length; index+=1) { 
            if(array[index]==OFF) { 
                sb.insert(bi, offInsert);
                bi += offLength;
            } else { 
                bi += 1;
            }
        }
        return sb.toString();
    }
    
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#isOn(int)
	 */
    public boolean isOn(int index) { return array[index] == ON; }
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#isOff(int)
	 */
    public boolean isOff(int index) { return array[index] == OFF; }
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#length()
	 */
    public int length() { return array.length; }

    /**
     * @param offset The offset up to which we should count (inclusive).
     * @return The number of "ON" bits with index <= offset.
     */
    public int countOnBits(int offset) { 
    	int count = 0; 
    	for(int i = 0; i <= offset; i++) { 
    		if(array[i] == ON) { 
    			count += 1;
    		}
    	}
    	return count;
    }
    
    /**
     * @param onBitIndex The number of the "on" bit which we want the index of.
     * @return The index (into the total bit-vector) of the specified "on" bit,
     * or -1 if onBitIndex < 0 or onBitIndex > # of ON bits.
     */
    public int reverseMap(int onBitIndex) {
    	int bitsSeen = 0;
    	int i = 0;
    	for(; i < array.length && bitsSeen <= onBitIndex; i++) { 
    		if(array[i] == ON) { bitsSeen += 1; }
    	}
    	if(i < array.length) { return i-1; }
    	return -1;
    }
    
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#countOnBits()
	 */
    public int countOnBits() { 
    	return countOnBits(array.length-1);
    }   
    
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#turnOnBit(int)
	 */
    public void turnOnBit(int index) { array[index] = ON; }
    /* (non-Javadoc)
	 * @see edu.mit.csail.cgs.utils.BitVector#turnOffBit(int)
	 */
    public void turnOffBit(int index) { array[index] = OFF; }
    
    public void bitwiseOR(ShortBitVector bv) { 
        if(array.length != bv.array.length) { throw new IllegalArgumentException(); }
        for(int i = 0; i < array.length; i++) { 
            if(array[i] == ON || bv.array[i] == ON) { 
                array[i] = ON; 
            } else { 
                array[i] = OFF;
            }
        }
    }
    
    public void bitwiseAND(ShortBitVector bv) { 
        if(array.length != bv.array.length) { throw new IllegalArgumentException(); }
        for(int i = 0; i < array.length; i++) { 
            if(array[i] == ON && bv.array[i] == ON) { 
                array[i] = ON; 
            } else { 
                array[i] = OFF;
            }
        }
    }
    
    public void flipLR() { 
        short[] newArray = new short[array.length];
        for(int i = 0; i < array.length; i++) { 
            newArray[array.length-i-1] = array[i];
        }
        array = newArray;
    }
    
    public boolean matches(String pattern) { 
        if(pattern.length() != array.length) { return false; }
        for(int i = 0; i < array.length; i++) { 
            char p = pattern.charAt(i);
            if(p != '?') { 
                if(p != '1' && p != '0') { return false; }
                if(p == '1' && array[i] == OFF) { return false; }
                if(p == '0' && array[i] == ON) { return false; }
            }
        }
        return true;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ShortBitVector)) { return false; }
        ShortBitVector bv = (ShortBitVector)o;
        if(array.length != bv.array.length) { return false; }
        for(int i = 0; i< array.length; i++) { 
            if(array[i] != bv.array[i]) { return false; }
        }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        for(int i = 0; i < array.length; i++) { 
            if(array[i] == ON) { 
                code += i; code *= 37;
            }
        }
        return code;
    }
    
    public String toString() { 
        return buildShortBitArrayString(array, 0, array.length);
    }
    
    public int compareTo(ShortBitVector bv) { 
        if(array.length < bv.array.length) { return -1; }
        if(array.length > bv.array.length) { return 1; }
        
        for(int i = 0; i < array.length; i++) { 
            if(array[i] < bv.array[i]) { return -1; }
            if(array[i] > bv.array[i]) { return 1; }
        }
        return 0;
    }
    
    /*
     * Static Helper Methods
     */
    
    public static Vector<ShortBitVector> createAllBitVectors(int length) { 
        Vector<ShortBitVector> vector = new Vector<ShortBitVector>();        
        ShortBitVector base = new ShortBitVector(length);
        addAllBitVectors(vector, base, 0);
        return vector;
    }
    
    private static void addAllBitVectors(Vector<ShortBitVector> vector, ShortBitVector base, int position) { 
        if(position < base.length()) { 
            base.turnOffBit(position);
            addAllBitVectors(vector, base, position+1);
            base.turnOnBit(position);
            addAllBitVectors(vector, base, position+1);
        } else { 
            vector.add(new ShortBitVector(base));
        }
    }

    private static int findHexChar(char c) { 
        for(int i = 0; i < hexChars.length; i++) { 
            if(hexChars[i] == c) { return i; }
        }
        return -1;
    }
    
    private static char translateShort2Hex(short[] array, int offset, int length) { 
        int hexIndex = 0;
        int placeSize = 1;
        for(int i = offset + length - 1; i >= offset; i--) { 
            if(array[i] == ON) { hexIndex += placeSize; }
            placeSize *= 2;
        }
        return hexChars[hexIndex];
    }
    
    private static void translateHex2Short(char hexChar, short[] array, int offset, int length) { 
        int value = findHexChar(hexChar);
        if(value == -1) { throw new IllegalArgumentException(); }
        int placeSize = 2;
        for(int i = offset+length-1; i>= offset; i--) {
            int rem = value % placeSize;
            if(rem != 0) { 
                array[i] = ON;
            } else { 
                array[i] = OFF;
            }
            value -= rem;
            placeSize *= 2;
        }
    }
    
    private static String buildShortBitArrayString(short[] array, int offset, int length) { 
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < length; i++) { 
            if(array[offset+i] == ON) { 
                sb.append("1");
            } else { 
                sb.append("0");
            }
        }
        return sb.toString();
    }
    
}
