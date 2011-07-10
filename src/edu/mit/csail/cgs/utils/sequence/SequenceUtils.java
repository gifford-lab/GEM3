package edu.mit.csail.cgs.utils.sequence;

import java.util.Collection;
import java.util.LinkedList;

/**
 * <code>SequenceUtils</code> provides a number of static methods for manipulating
 * DNA sequences stores as strings or char[].
 *
 * @author <a href="mailto:arolfe@mit.edu">Alex Rolfe</a>
 */
public class SequenceUtils {

    
    /**
     * <code>complement</code> returns the complement of a nucleotide in the 2-bit representation.
     */
    public static int complement(int i) { 
    	switch(i) { 
    	case 0: return 1;
    	case 1: return 0;
    	case 2: return 3; 
    	case 3: return 2;
    	default: return -1;
    	}
    }
    /**
     * <code>complement</code> returns the complement of a nucleotide in the character representation
     * (A,C,T,G,a,c,t,g).
     */
    public static char complementChar(char c) { 
        if (trans == null) {
            trans = new char['z'];
            trans['A'] = 'T';
            trans['C'] = 'G';
            trans['T'] = 'A';
            trans['G'] = 'C';
            trans['a'] = 'T';
            trans['c'] = 'G';
            trans['t'] = 'A';
            trans['g'] = 'C';
            trans['n'] = 'N';
            trans['N'] = 'N';
        }
        return trans[c];
    }
    


    private static char[] trans;
    /**
     * <code>reverseComplement</code> mutates in the input array of characters
     * (A,C,T,G,a,c,t,g) to be the reverse complement.
     */
    public static void reverseComplement(char[] array) {
        if (trans == null) {
            trans = new char['z'];
            trans['A'] = 'T';
            trans['C'] = 'G';
            trans['T'] = 'A';
            trans['G'] = 'C';
            trans['a'] = 't';
            trans['c'] = 'g';
            trans['t'] = 'a';
            trans['g'] = 'c';
            trans['N'] = 'N';
            trans['n'] = 'n';
            trans['X'] = 'X';
            trans['x'] = 'x';

        }
        int i;
        int end = array.length - 1;        
        for (i = 0; i <= array.length / 2 && i < array.length; i++) {
            try {
                char first = array[i];
                array[i] = trans[array[end - i]];
                array[end-i] = trans[first];
            } catch (ArrayIndexOutOfBoundsException ex) {
                ex.printStackTrace();
                System.err.println("i=" + i);
                System.err.println("first = " + array[i]);
                System.err.println("other = " + array[end-i]);
                System.err.println("trans = " + trans[array[i]] + " and " + trans[array[end-i]]);
            }
        }
    }

    public static String reverseComplement(String str) { 
		StringBuilder sb = new StringBuilder();
		for(int i = str.length()-1; i>= 0; i--) { 
			sb.append(complementChar(str.charAt(i)));
		}
		return sb.toString();
	}
    
    public static byte[] reverseComplement(byte[] bases) {
		byte[] rc = new byte[bases.length];
		int j=0;
		for(int i = bases.length-1; i>= 0; i--) { 
			rc[j]=(byte)complementChar((char)bases[i]);
			j++;
		}
		return rc;
	}

    /** converts from 2-bit representation to character representation
     */
    public static char int2Char(int i) { 
        switch(i) { 
        case 0: return 'A';
        case 1: return 'T';
        case 2: return 'G';
        case 3: return 'C';
        }
        return 'n';
    }
    
    /** converts from character representation to 2-bit representation
     */
    public static int char2Int(char c) { 
        switch(c) { 
        case 'a':
        case 'A': 
            return 0; 
        case 't':
        case 'T':
            return 1;
        case 'g':
        case 'G':
            return 2;
        case 'c':
        case 'C':
            return 3;
        }
        return -1;
    }

    /** maps a DNA sequence (ie [actgACTG]*) to 
     * a Long representing that sequence in 2-bit per
     * base representation
     *
     *
     * we should time this and compare performance to the version
     * where you have 
     * int charToInt[256]
     * charToInt['a'] = 0
     * etc
     * and do lookups from that...
     */
    public static long StringToLong(String a) {
        return StringToLong(a,0,a.length());
    }

    /** convert k characters of String to a long, starting at position offset
     * 
     * Throws a StringIndexOutOfBoundsException (or something) if there aren't
     * enough characters in the string
     *
     * Apparently Tim and Alex did this conversion in opposite order- Tim's had
     * string[0] as the lowest bits in the long and Alex had them as the highest.  
     * I like my way better so that
     *  1) when you write out the long and the string, the letters are in the same order
     *  2) the strings and the longs sort the same order
     *
     */
    public static long StringToLong(String a, int offset, int k) {
        long sum = (long)0;
        for(int i = 0; i < k; i++) { 
            int val = char2Int(a.charAt(i + offset));
            sum = (sum << 2) + val;
        }
        return sum;
    }

    public static String LongToString(Long l, int length) {
        char[] output = new char[length];
        while (length-- > 0) {
            output[length] = int2Char((int)(l & 3));
            l >>= 2;
        }
        return new String(output);
    }

	public static Collection<String> creatKMers(String original, int k) { 
		LinkedList<String> kmers = new LinkedList<String>();
		for(int i = 0; i <= original.length() - k; i++) { 
			kmers.addLast(original.substring(i, i + k));
		}
		return kmers;
	}    
	
	private static cern.jet.random.engine.RandomEngine randomEngine;
	public static String generateRandomBases(int k){
		if (randomEngine==null)
			randomEngine = new cern.jet.random.engine.MersenneTwister();
		StringBuffer sb=new StringBuffer();
		for (int i=0;i<k;i++){
			int num = Math.abs(randomEngine.nextInt()) % 4;
			sb.append(int2Char(num));
		}
		return sb.toString();
	}
	public static String generateRandomString(int k){
		if (randomEngine==null)
			randomEngine = new cern.jet.random.engine.MersenneTwister();
		char[] chars = new char[k];
		for (int i=0;i<k;i++){
			chars[i]=(char)(Math.abs(randomEngine.nextInt()) % 26 + 'A');
		}
		return new String(chars);
	}

}
