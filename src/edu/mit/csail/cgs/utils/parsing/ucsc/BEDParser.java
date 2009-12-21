package edu.mit.csail.cgs.utils.parsing.ucsc;

import java.util.*;
import java.util.regex.*;
import java.awt.Color;
import java.io.*;

/**
 * @author tdanford
 * 
 * BED Format: 
 * http://genome.ucsc.edu/FAQ/FAQformat#format1
 * 
 * The first three required BED fields are:

   1. chrom - 
   The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
   2. chromStart - The starting position of the feature in the chromosome or scaffold. 
   The first base in a chromosome is numbered 0.
   3. chromEnd - The ending position of the feature in the chromosome or scaffold. 
   The chromEnd base is not included in the display of the feature. For example, 
   the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, 
   and span the bases numbered 0-99. 

The 9 additional optional BED fields are:

   4. name - Defines the name of the BED line. This label is displayed to the left of 
   the BED line in the Genome Browser window when the track is open to full display 
   mode or directly to the left of the item in pack mode.
   5. score - A score between 0 and 1000. If the track line useScore attribute is set 
   to 1 for this annotation data set, the score value will determine the level of gray 
   in which this feature is displayed (higher numbers = darker gray).
   6. strand - Defines the strand - either '+' or '-'.
   7. thickStart - The starting position at which the feature is drawn thickly 
   (for example, the start codon in gene displays).
   8. thickEnd - The ending position at which the feature is drawn thickly (for 
   example, the stop codon in gene displays).
   9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line 
   itemRgb attribute is set to "On", this RBG value will determine the display color 
   of the data contained in this BED line. NOTE: It is recommended that a simple 
   color scheme (eight colors or less) be used with this attribute to avoid overwhelming 
   the color resources of the Genome Browser and your Internet browser.
  10. blockCount - The number of blocks (exons) in the BED line.
  11. blockSizes - A comma-separated list of the block sizes. The number of items in 
  this list should correspond to blockCount.
  12. blockStarts - A comma-separated list of block starts. All of the blockStart 
  positions should be calculated relative to chromStart. The number of items in 
  this list should correspond to blockCount. 
  
 */
public class BEDParser implements Iterator<BEDLine> {
	
	private BEDLine nextLine;
	private BufferedReader br;
	
	public BEDParser(File f) throws IOException { 
		br = new BufferedReader(new FileReader(f));
		findNextLine();
	}
	
	private void findNextLine() throws IOException { 
		String line = br.readLine();
		if(line != null) { 
			nextLine = new BEDLine(line);
		} else { 
			nextLine = null;
			br.close();
		}
	}
	
	public boolean hasNext() { 
		return nextLine != null;
	}
	
	public BEDLine next() { 
		BEDLine ret = nextLine;
		try {
			findNextLine();
		} catch (IOException e) {
			e.printStackTrace();
			nextLine = null;
			try {
				br.close();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
		}
		return ret;
	}
	
	public void remove() { 
		throw new UnsupportedOperationException();
	}
	
}
