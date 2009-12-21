/*
 * Author: tdanford
 * Date: Aug 24, 2008
 */
package edu.mit.csail.cgs.datasets.alignments;

import java.util.HashMap;
import java.util.Vector;

import edu.mit.csail.cgs.utils.ArrayUtils;
import edu.mit.csail.cgs.utils.Predicate;

/**
 * A GappedAlignmentString is part of a two-way (or multiple) alignment. 
 * Each object of GappedAlignemntString represents a single "gapped alignment," 
 * which is a sequence that is interspersed with gap symbols ('-') that are 
 * used to indicate where other aligned sequences contain insertions relative
 * to this sequence.  
 * 
 * Gapped alignment strings are usually written "on top" of each other, in an alignment. 
 * For example: 
 * AATGAC---CAG
 * AA--ACGGGGAG
 * 
 * This is two (different) gapped alignmemnt strings.  Each gapped alignment string 
 * implicitly uses two different coordinate spaces: "gapped" and "ungapped" coordinates.
 * 
 * Take the second alignment string: "AA--ACGGGGAG".  
 * 
 * The "gapped coordinates" index characters from this string. So, in ungapped coordinates,
 * 0: A
 * 1: A
 * 2: - 
 * 3: -
 * 4: A 
 * and so on.  
 * 
 * The "ungapped coordinates," by contrast, index only the original string without gaps: 
 * 0: A
 * 1: A
 * 2: A
 * 3: C
 * 4: G
 * etc.  
 * 
 * GappedAlignemntString provides a way of converting between these two coordinate 
 * systems, for a particular gapped alignment string.  It is used extensively as part of the 
 * MultipleAlignment class, as well.   
 * 
 * @author tdanford
 */
public class GappedAlignmentString {
	
	private char[] gapped;
	private int ungappedLength;  // assertion: ungappedLength <= gapped.length
	
	/**
	 * The gap positions
	 */
	private int[] gapPositions;
	
	/**
	 * Mapping of the sequence's position (in ungapped form) to its position (in gapped form)
	 */
	private int[] ungappedMap;
	
	/**
	 * Mapping of the sequence's position (in gapped form) to its position (in ungapped form)</br>
	 * Note: if the position (in gapped form) corresponds to a gap (dash, -), then it is assigned
	 * the position (in ungapped form) which is closer to its right.
	 */
	private int[] gappedMap;

	public GappedAlignmentString(char[] array) { 
		gapped = array.clone();
		ungappedLength = 0;
		for(int i = 0; i < gapped.length; i++) {
			if(gapped[i] != '-') { 
				ungappedLength += 1;
			}
		}

		//gapPositions = determineGapPositions();
		gapPositions = new int[gapped.length-ungappedLength];
		
		gappedMap = new int[gapped.length];
		ungappedMap = new int[ungappedLength];
		for(int i = 0, j = 0, k = 0; i < gapped.length; i++) { 
			gappedMap[i] = j;
			if(gapped[i] != '-') {
				ungappedMap[j++] = i;
			} else { 
				gapPositions[k++] = i;
			}
		}
		
	}
	
	/**
	 * Returns a GappedAlignmentString built from a substring of the original 
	 * gapped alignment.  The substring is taken from the range [start, end) of the 
	 * original string.  
	 * @param start Starting coordinate (inclusive)
	 * @param end Ending coordinate (exclusive)
	 * @return
	 */
	public GappedAlignmentString substring(int start, int end) { 
		char[] array = new char[end-start];
		for(int i = start; i < end; i++) { 
			array[i-start] = gapped[i];
		}
		return new GappedAlignmentString(array);
	}
	
	/**
	 * 
	 * @return The gapped length for the sequence.
	 */
	public int gappedLength() { return gapped.length; }
	
	/**
	 * 
	 * @return The ungapped length for the sequence.
	 */
	public int ungappedLength() { return ungappedLength; }
	
	/**
	 * 
	 * @param i The gapped coordinate to find its char
	 * @return The char that this position corresponds to
	 */
	public char gappedChar(int i) { return gapped[i]; }

	public String gappedSequence(int j1, int j2) {
		char[] seq = new char[j2-j1+1];
		for(int i = j1, k = 0; i <= j2; i++, k++) { 
			seq[k] = gapped[i];
		}
		return new String(seq);
	}

	/**
	 * 
	 * @param i The ungapped coordinate to find its char
	 * @return The char that this position corresponds to
	 */
	public char ungappedChar(int i) { return gapped[ungappedMap[i]]; }
	
	/**
	 * Converts from gapped to ungapped coordinates.  
	 * Therefore, in the gapped alignment string A--AGCA, 
	 * gappedToUngapped(0) = 0, and
	 * gappedToUngapped(1) = 0, and 
	 * gappedToUngapped(3) = 1
	 * 
	 * @param i  The gapped coordinate to convert.
	 * @return
	 */
	public int gappedToUngapped(int i) { return gappedMap[i]; }

	/**
	 * Converts from ungapped to gapped coordinates.  
	 * Therefore, in the gapped alignment string A--AGCA, 
	 * ungappedToGapped(0) = 0
	 * ungappedToGapped(1) = 3
	 * ungappedToGapped(2) = 4
	 * 
	 * @param i  The ungapped coordinate to convert.
	 * @return
	 */
	public int ungappedToGapped(int i) { return ungappedMap[i]; }
	
	/**
	 * @param gi The (gapped) coordinate.
	 * @return <tt>true</tt> if there is a gap in this position.</br>
	 * <tt>false</tt> otherwise 
	 */
	public boolean isGap(int gi) { return gapped[gi] == '-'; }
	
	/**
	 * 
	 * @return The ungapped sequence as a <tt>String</tt>.
	 */
	public String ungappedString() { 
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < ungappedLength; i++) { 
			sb.append(ungappedChar(i));
		}
		return sb.toString();
	}

	/**
	 * 
	 * @return The gapped sequence as a <tt>String</tt>.
	 */
	public String gappedString() {
		return new String(gapped);
	}
	
	/**
	 * 
	 * @return a mapping of gapped position to its corresponding ungapped position
	 * </br>Contains the same info as <tt>gappedMap</tt>, however as a <tt>HashMap</tt>
	 * variable it allows O(1) access.
	 * @deprecated Prefer using <tt>gappedToUngapped</tt> method
	 * @see GappedAlignmentString#gappedToUngapped(int)
	 * @see HashMap
	 */
	@Deprecated
	public HashMap<Integer, Integer> getGapped2UngappedMap()
	{
		int N = gappedMap.length;
		
		// making map allocation more efficient
		int initCapacity  = (int)Math.ceil(1.5*N);
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>(initCapacity);
		
		for(int i = 0; i < N; i++)
			map.put(i, gappedMap[i]);
		
		return map;
		
	}//end of getGapped2UngappedMap method
	
	/**
	 * @return a mapping of gapped position to its corresponding ungapped position
	 * </br>Contains the same info as <tt>gappedMap</tt>, however as a <tt>HashMap</tt>
	 * variable it allows O(1) access.
	 * @deprecated Prefer using <tt>ungappedToGapped</tt> method
	 * @see GappedAlignmentString#ungappedToGapped(int)
	 * @see HashMap
	 */
	@Deprecated
	public HashMap<Integer, Integer> getUngapped2GappedMap()
	{
		int N = ungappedMap.length;
		
		// making map allocation more efficient
		int initCapacity  = (int)Math.ceil(1.5*N);
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>(initCapacity);
		
		for(int i = 0; i < N; i++)
			map.put(i, ungappedMap[i]);
		
		return map;
		
	}//end of getUngapped2GappedMap method
	
	/**
	 * 
	 * @return The gap positions in this sequence
	 */
	public Integer[] determineGapPositions()
	{
		Integer[] indices = ArrayUtils.range(0, gapped.length);
		return ArrayUtils.mask(indices, new Predicate<Integer>() {
			public boolean accepts(Integer v) {
				return isGap(v);
			} 
		});
		/*
		Vector<Integer> gapPositionsVector = new Vector<Integer>();
		
		int N = gapped.length;
		for(int i = 0; i < N; i++)
		{
			if( isGap(i) )
				gapPositionsVector.add(i);
		}
		
		Integer[] gapPositions = gapPositionsVector.toArray(new Integer[gapPositionsVector.size()]);
		
		return gapPositions;
		*/
	}//end of determineGapPositions method
	
	public int[] getGapPositions()
	{
		return gapPositions;
	}//end of getGapPositions method

}
