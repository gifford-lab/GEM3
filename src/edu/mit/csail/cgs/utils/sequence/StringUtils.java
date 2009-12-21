package edu.mit.csail.cgs.utils.sequence;

import java.util.*;

/**
 * @author tdanford
 * StringUtils: a set of utility methods for dealing with biological sequences and strings.
 * 
 * Back in the day, we used to use the old BioJava DNAUtils class and its correlates to do 
 * all that.  Since we've drifted away from that framework, we have needed to roll our own
 * implementations of the same framework.  But really, we should all be using the same one,
 * so let's put it here and use *this* one.
 *
 * @deprecated Use SequenceUtils instead since there's a cgs.utils.strings.StringUtils too and
 *  that's confusing
 */
public abstract class StringUtils { 

	public static String revcomp(String str) { 
		StringBuilder sb = new StringBuilder();
		for(int i = str.length()-1; i>= 0; i--) { 
			sb.append(complement(str.charAt(i)));
		}
		return sb.toString();
	}
	
	public static char complement(char c) { 
		switch(c) { 
		case 'A': return 'T';
		case 'T': return 'A';
		case 'G': return 'C';
		case 'C': return 'G';
		case 'a': return 't';
		case 't': return 'a';
		case 'g': return 'c';
		case 'c': return 'g';
		default: return c;
		}
	}
	
	public static Collection<String> creatKMers(String original, int k) { 
		LinkedList<String> kmers = new LinkedList<String>();
		for(int i = 0; i <= original.length() - k; i++) { 
			kmers.addLast(original.substring(i, i + k));
		}
		return kmers;
	}
}
