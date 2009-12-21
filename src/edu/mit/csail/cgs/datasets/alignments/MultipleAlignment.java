/*
 * Author: tdanford
 * Date: Aug 24, 2008
 */
package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.util.HashMap;


import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.utils.strings.StringUtils;

/**
 * Ties several GappedAlignmentStrings into a single "multiple alignment," which is
 * guaranteed to have the same gapped length across all strings.  
 * 
 * Gives column-wise access to the multiple alignment, as well as other utility functions.  
 * 
 * @author tdanford
 * @see GappedAlignmentString
 *
 */
public class MultipleAlignment {

	/**
	 * A map containing the names of the sequences and the sequences themselves (as 
	 * <tt>GappedAlignmentStrings</tt>).
	 */
	private Map<String,GappedAlignmentString> gappedAlignments;
	
	/**
	 * The names of the sequences as a <tt>Vector</tt> of <tt>Strings</tt>
	 */
	private Vector<String> ordering;
	
	/**
	 * The gapped length of all sequences. </br>
	 * Same for all sequences.</br>
	 * (Dashes (-) are included in the counting)
	 */
	private int gappedLength;
	
	/**
	 * Creates a "blank" MultipleAlignment.  
	 *
	 */
	public MultipleAlignment() { 
		gappedAlignments = new HashMap<String,GappedAlignmentString>();
		ordering = new Vector<String>();
		gappedLength = 0;
	}
	
	/**
	 * Creates a MultipleAlignment from the given Alignment object.  Alignment is 
	 * roughly analogous to MultipleAlignment, but it's an interface and it corresponds 
	 * to a multiple alignment that has been read from a data store somewhere 
	 * (presumably a database, but perhaps also a file or other source).  Creating a 
	 * MultipleAlignment object "around" it basically makes it useable -- MultipleAlignment
	 * will read the Alignment object and its corresponding AlignBlock pieces, and 
	 * builds data structures so that the underlying multiple alignment can be 
	 * accessed easily.   
	 * 
	 * @param a
	 */
	public MultipleAlignment(Alignment a) { 
		this();
		Collection<AlignBlock> blocks = a.getAlignBlocks();
		for(AlignBlock block : blocks) { 
			GappedAlignmentString gaps = new GappedAlignmentString(block.getBitString());
			Genome g = block.getGenome();
			String chrom = block.getChrom();
			String name = String.format("%s_%s", g.getVersion(), chrom);
			addGappedAlignment(name, gaps);
		}
	}
	
	private MultipleAlignment(Vector<String> ord, String[] seqs) { 
		this();
		for(int i = 0; i < ord.size(); i++) { 
			GappedAlignmentString gap = new GappedAlignmentString(seqs[i].toCharArray());
			addGappedAlignment(ord.get(i), gap);
		}
	}
	
	/**
	 * Adds the sequence with name <tt>name</tt> and sequence <tt>gas</tt> to the alignment.
	 * 
	 * @param name The name of the sequence to be added 
	 * @param gas The sequence as a <tt>GappedAlignmentString</tt>
	 * @see		GappedAlignmentString
	 */
	public void addGappedAlignment(String name, GappedAlignmentString gas) { 
		if(ordering.isEmpty()) { 
			gappedLength = gas.gappedLength();
		} else { 
			if(gas.gappedLength() != gappedLength) { 
				throw new IllegalArgumentException(String.format("Gapped Length %d doesn't match.", 
						gas.gappedLength()));
			}
			if(ordering.contains(name)) { 
				throw new IllegalArgumentException(String.format("Duplicate name: %s", name)); 
			}
		}
		ordering.add(name);
		gappedAlignments.put(name, gas);
	}
	
	/**
	 * @return The gapped length for all sequences. </br>
	 * Note: Same for all
	 */
	public int gappedLength() { return gappedLength; }
	
	/**
	 * 
	 * @return The number of species in the multiple alignment. 
	 */
	public int numSpecies() { return ordering.size(); }
	
	/**
	 * 
	 * @return The species names as <tt>Strings</tt>
	 */
	public String[] species() { return ordering.toArray(new String[ordering.size()]); }
	
	/**
	 * 
	 * @param str Name of the sequence whose ungapped length is queried 
	 * @return The ungapped length of this sequence
	 */
	public int ungappedLength(String str) { return gappedAlignments.get(str).ungappedLength(); }
	
	/**
	 * 
	 * @param spec Name of the sequence whose ungapped sequence is queried
	 * @return The ungapped sequence as a <tt>String</tt>
	 */
	public String ungappedString(String spec) { return gappedAlignments.get(spec).ungappedString(); }
	
	/**
	 * 
	 * @param spec Name of the sequence whose gapped sequence is queried
	 * @return The gapped sequence as a <tt>String</tt>
	 */
	public String gappedString(String spec) { return gappedAlignments.get(spec).gappedString(); }
	
	/**
	 * 
	 * @param spec Name of the sequence whose sequence is queried
	 * @return The sequence as a <tt>GappedAlignmentString</tt> object
	 * @see GappedAlignmentString
	 */
	public GappedAlignmentString getGappedAlignment(String spec) { 
		return gappedAlignments.get(spec);
	}
	
	/**
	 * Checks if there is a gap in position <tt>gi</tt> in at least one of the 
	 * aligned sequences
	 * @param gi Position to be queried for the presence of a gap
	 * @return <tt>true</tt> if there is no gap in any of the aligned sequences</br>
	 * <tt>false</tt> if there is at least one. 
	 */
	public boolean isPresent(int gi) { 
		for(int i = 0; i < ordering.size(); i++) { 
			if(gappedAlignments.get(ordering.get(i)).isGap(gi)) { 
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Returns whether all columns are "present" (i.e. non-gapped) in a set of indices, 
	 * not just a single column. 
	 * 
	 * @param g1
	 * @param g2
	 * @return
	 */
	public boolean isPresent(int g1, int g2) { 
		for(int i = g1; i <= g2; i++) { 
			if(!isPresent(i)) { 
				return false;
			}
		}
		return true;
	}

	private int maxTagLength() { 
		int tl = 0;
		for(String tag : gappedAlignments.keySet()) { 
			tl = Math.max(tl, tag.length());
		}
		return tl;
	}
	
	/**
	 * returns the multiple alignment as a string in a fashion 
	 * identical to that in a file of ALN/CLUSTALW format (except from
	 * the fact that info about residues' match is missing. Aka, the 3rd
	 * line of a file of ALN/CLUSTALW ('*', '.', ':' chars) is missing).
	 * 
	 * (Note: CLUSTALW doesn't necessarily have the match characters... -Tim)
	 */
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		int tl = maxTagLength();
		int i = 0, last = gappedAlignments.keySet().size()-1;
		for(String tag : gappedAlignments.keySet()) { 
			sb.append(StringUtils.padString(tag, tl));
			sb.append("  ");
			sb.append(gappedAlignments.get(tag).gappedString());
			if(i < last) { sb.append("\n"); }
			i += 1;
		}
		return sb.toString();
	}
	
	/**
	 * @param j  The ungapped index into the multiple alignment.
	 * @return An array of characters, corresponding to a particular column from the multiple
	 * alignment.
	 */
	public char[] alignedChars(int j) { 
		char[] array = new char[ordering.size()];
		for(int i = 0; i < ordering.size(); i++) {
			String species = ordering.get(i);
			array[i] = gappedAlignments.get(species).gappedChar(j);
		}
		return array;
	}
	
	/**
	 * Returns a section of the multiple alignment in the region between the gapped indices
	 * [j1, j2] (notice that the coordinates are inclusive).  
	 * 
	 * @param j1 The start coordinate (inclusive)
	 * @param j2 The end-coordinate (inclusive)
	 * @return An array of strings, one for each sequence in the multiple alignment. 
	 */
	public String[] alignedSeqs(int j1, int j2) {
		if(j2 < j1) { 
			throw new IllegalArgumentException(String.format("[%d, %d] are illegal coordinates " +
					"for MultipleAlignment.alignedSeqs()", j1, j2));
		}
		
		String[] array = new String[ordering.size()];
		for(int i = 0; i < ordering.size(); i++) { 
			String species = ordering.get(i);
			array[i] = gappedAlignments.get(species).gappedSequence(j1, j2);
		}
		return array;
	}
	
	/**
	 * Returns a new MultipleAlignment object over the same set of sequences, but which 
	 * contains *only* the gapped alignments between the given indices [j1, j2].  
	 * 
	 * @param j1
	 * @param j2
	 * @return
	 */
	public MultipleAlignment subAlignment(int j1, int j2) { 
		return new MultipleAlignment(ordering, alignedSeqs(j1, j2));
	}
	
	public int mapUngapped(String s1, String s2, int ugidx) { 
		return mapUngapped(ordering.indexOf(s1), ordering.indexOf(s2), ugidx);
	}
	
	
	/**
	 * Valid values for <tt>s1</tt> and <tt>s2</tt> are [0, ordering.size()-1]
	 * @param s1 The serial number of the sequence (<tt>s1</tt>) whose coordinate we wish to map against 
	 * that of (<tt>s2</tt>)
	 * @param s2 The other sequence
	 * @param ugidx The ungapped coordinate of (<tt>s1</tt>)
	 * @return The ungapped coordinate of (<tt>s2</tt>) corresponding to <tt>ugidx</tt>
	 * @see  GappedAlignmentString
	 */
	public int mapUngapped(int s1, int s2, int ugidx) throws IndexOutOfBoundsException 
	{ 
		if( ( s1 < 0 || s1 >= ordering.size() ) || ( s2 < 0 || s2 >= ordering.size() ))
			throw new IndexOutOfBoundsException("Valid values for s1 and s2 are [0, ordering.size()-1]");
		
		GappedAlignmentString str1 = gappedAlignments.get(ordering.get(s1));
		GappedAlignmentString str2 = gappedAlignments.get(ordering.get(s2));
		int gidx = str1.ungappedToGapped(ugidx);
		int ugidx2 = str2.gappedToUngapped(gidx);
		if(str2.isGap(gidx)) { return -1; }
		return ugidx2;
	}

	/*
	 * TODO: Isn't all the information assembled by the next three functions already 
	 * available through the GappedAlignmentString class? 
	 */

	/**
	 * This function returns a map whose entries hold the coordinates of the sequence
	 * as keys and the coordinates of the target sequences as values.</br>
	 * If a coordinate in sequence corresponds to a gap in another sequence, then -1 is
	 * assigned to that position.</br>
	 * For example suppose we have the sequences:</br>
	 * seq1: AAT-GCT-TCG</br>
	 * seq2: G--C-CTCT-C</br>
	 * seq3: A-T--GT-AG-</br>
	 * Then, if we wish to find the mapping of seq1's coordinates to seq2's coordinates, this will
	 * give us: 0 -> 0, 1 -> -1, 2 -> -1, 3 -> -1, 4 -> 2, 5 -> 3, 6 -> 5, 7 -> -1, 8 -> 6.</br>
	 * Also, if we wish to find the mapping of seq1's coordinates to seq2's and seq3's coordinates, this will
	 * give us: 0 -> {0, 0}, 1 -> {-1, -1}, 2 -> {-1, 1}, 3 -> {-1, -1}, 4 -> {2, 2}, 5 -> {3, 3}, 6 -> {5, 4}, 7 -> {-1, 5}, 8 -> {6, -1}.
	 * @param <T> This could be either the serial number of the sequence (valid range: [0, sequences.length-1])
	 * or the sequence names 
	 * @param <E>
	 * @param seqID The sequence ID (either serial number or name) which we wish to map</br>
	 * Note that the valid range of serial numbers is [0, sequences.length-1]
	 * @param targetSeqIDs The (target) sequences IDs (either serial numbers or names)
	 * @return A map whose keys are the sequence coordinates and values the corresponding coordinates of the target
	 * sequences 
	 * @throws IllegalArgumentException
	 * @see GappedAlignmentString
	 */
	public <T extends Object>   Map mapSeqCoords2SeqsCoords(T generic_seqID, Vector<T> generic_targetSeqIDs)
	throws IllegalArgumentException, IndexOutOfBoundsException
	{
		GappedAlignmentString gas_seq;
		GappedAlignmentString[] gas_targetSeqs;
		
		String[] speciesNames = species();
		
		// Check if generic_seqID is mistakenly included in the generic_targetSeqIDs set 
		if( generic_targetSeqIDs.contains(generic_seqID) )
			generic_targetSeqIDs.remove(generic_seqID);
			
		if( generic_targetSeqIDs.size() + 1 >  speciesNames.length )
			throw new IndexOutOfBoundsException("You entered more sequences than are in the file.");
	
		Map map = new HashMap();
		// Check whether sequence IDs are inputed as serial numbers or as names
		String seqsClassName = generic_seqID.getClass().getSimpleName();
		
		//if( generic_seqID.getClass() instanceof java.lang.Integer )
			//System.out.println("XAXA");
		
		if( seqsClassName.equals("Integer") )
		{
			Integer seqID = (Integer)generic_seqID;
			Integer[] seqIDs = generic_targetSeqIDs.toArray(new Integer[generic_targetSeqIDs.size()]);
		
			if( (StatUtil.findMax(seqIDs).getFirst() > speciesNames.length -1) || (StatUtil.findMin(seqIDs).getFirst() < 0) )
				throw new IndexOutOfBoundsException("The sequence IDs (aka, serial numbers have to lie inside" + 
                                                    " the range 0 to speciesNames.length -1");

			gas_seq = getGappedAlignment(speciesNames[seqID]);
			gas_targetSeqs = new GappedAlignmentString[seqIDs.length];
			for(int i = 0; i < seqIDs.length; i++)
				gas_targetSeqs[i] = getGappedAlignment(speciesNames[seqIDs[i]]);
		}
		else if( seqsClassName.equals("String") )
		{
			String seqID = (String)generic_seqID;
			String[] seqIDs = generic_targetSeqIDs.toArray(new String[generic_targetSeqIDs.size()]);
	
			gas_seq = getGappedAlignment(seqID);
			gas_targetSeqs = new GappedAlignmentString[seqIDs.length];
			for(int i = 0; i < seqIDs.length; i++)
				gas_targetSeqs[i] = getGappedAlignment(seqIDs[i]);
		}
		else
		{
			throw new IllegalArgumentException("Sequence IDs should be either of type Integer or String");
		}
		
		// Pairwise Alignment
		if( gas_targetSeqs.length == 1 )
		{
			map = doMapSeq2Seq(gas_seq, gas_targetSeqs);
		}
		// Multiple Alignment
		else
		{
			map = doMapSeq2Seqs(gas_seq, gas_targetSeqs);
		}
		
		return map;
		
	}//end of mapSeqCoords2SeqsCoords method
	
	/**
	 * 
	 * @param gas_seq The sequence which we wish to map it against other sequence(s)
	 * @param gas_targetSeqs The target sequence
	 */
	private HashMap<Integer, Integer> doMapSeq2Seq(GappedAlignmentString gas_seq, GappedAlignmentString[] gas_targetSeqs)
	{
		int N = gas_seq.gappedLength();
		// making map allocation more efficient
		int initCapacity  = (int)Math.ceil(1.5*N);
		HashMap<Integer, Integer>  map = new HashMap<Integer, Integer>(initCapacity); 
		
		//When a gap position corresponds to a gap remove it from seq
		HashMap<Integer, Integer> gapped2UngappedMap_seq = gas_seq.getGapped2UngappedMap();
		Integer[] gapPositions_seq = gas_seq.determineGapPositions();
		for(Integer e:gapPositions_seq)
			gapped2UngappedMap_seq.remove(e);
		
		//When a gap position corresponds to a gap, put -1 as its value
		HashMap<Integer, Integer> gapped2UngappedMap_targetSeq = gas_targetSeqs[0].getGapped2UngappedMap();
		Integer[] gapPositions_targetSeq = gas_targetSeqs[0].determineGapPositions();
		for(Integer e:gapPositions_targetSeq)
			gapped2UngappedMap_targetSeq.put(e, -1);
		
		
		Integer[] gappedPositionsNoGaps_seq = gapped2UngappedMap_seq.keySet().toArray(new Integer[gapped2UngappedMap_seq.size()]);
		for(Integer gappedPosNoGaps : gappedPositionsNoGaps_seq)
		{
			Integer ungappedPos = gapped2UngappedMap_targetSeq.get(gappedPosNoGaps);
			if(ungappedPos.equals(null)){ ungappedPos = -1;}
			map.put(gapped2UngappedMap_seq.get(gappedPosNoGaps), ungappedPos);
		}
		
		return map;
	}//end of doMapSeq2Seq method
	
	/**
	 * 
	 * @param <E> Integer in case of pairwise alignment, ArrayList otherwise
	 * @param gas_seq The sequence which we wish to map it against other sequence(s)
	 * @param gas_targetSeqs The target sequences
	 * @return
	 */
	private   HashMap<Integer, Integer[]> doMapSeq2Seqs(GappedAlignmentString gas_seq, GappedAlignmentString[] gas_targetSeqs)
	{
		int N = gas_seq.gappedLength();
		// making map allocation more efficient
		int initCapacity  = (int)Math.ceil(1.5*N);
		HashMap<Integer, Integer[]>  map = new HashMap<Integer, Integer[]>(initCapacity); 
		
		//When a gap position corresponds to a gap remove it from seq
		HashMap<Integer, Integer> gapped2UngappedMap_seq = gas_seq.getGapped2UngappedMap();
		Integer[] gapPositions_seq = gas_seq.determineGapPositions();
		for(Integer e:gapPositions_seq)
			gapped2UngappedMap_seq.remove(e);
		
		//When a gap position corresponds to a gap, put -1 as its value
		HashMap<Integer, Integer>[] gapped2UngappedMap_targetSeqs = new HashMap[gas_targetSeqs.length];
		for(int i = 0; i < gas_targetSeqs.length; i++)
			gapped2UngappedMap_targetSeqs[i] = new HashMap<Integer, Integer>();
		
		int i = 0;
		for(GappedAlignmentString gas_targetSeq: gas_targetSeqs)
		{
		HashMap<Integer, Integer> gapped2UngappedMap_targetSeq = gas_targetSeq.getGapped2UngappedMap();
			Integer[] gapPositions_targetSeq = gas_targetSeq.determineGapPositions();
			for(Integer e:gapPositions_targetSeq)
				gapped2UngappedMap_targetSeq.put(e, -1);
			
			gapped2UngappedMap_targetSeqs[i++] = gapped2UngappedMap_targetSeq;
		}
		
		Integer[] gappedPositionsNoGaps_seq = gapped2UngappedMap_seq.keySet().toArray(new Integer[gapped2UngappedMap_seq.size()]);
		for(Integer gappedPosNoGaps : gappedPositionsNoGaps_seq)
		{
			Integer[] ungappedPos_seqs = new Integer[gas_targetSeqs.length];
			
			for(i = 0; i < gas_targetSeqs.length; i++)
			{
				ungappedPos_seqs[i] = gapped2UngappedMap_targetSeqs[i].get(gappedPosNoGaps);
				if(ungappedPos_seqs[i].equals(null)){ ungappedPos_seqs[i] = -1;}
			}
			
			map.put(gapped2UngappedMap_seq.get(gappedPosNoGaps), ungappedPos_seqs);
		}
		
		return map;
	}//end of doMapSeq2Seqs method 
}
