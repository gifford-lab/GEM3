package edu.mit.csail.cgs.utils.io.parsing.ucsc;

import java.awt.Color;
import java.util.regex.*;

public class BEDLine {
	
	private String chrom;
	private int chromStart, chromEnd;
	private String name;
	private Double score;
	private Character strand;
	private Integer thickStart, thickEnd;
	private Color itemRgb;
	private int blockCount;
	private int[] blockSizes;
	private int[] blockStarts;
	
	public String getChrom() { return chrom; }
	public int getChromStart() { return chromStart; }
	public int getChromEnd() { return chromEnd; }
	public String getName() { return name; }
	public Double getScore() { return score; }
	public Character getStrand() { return strand; }
	public Integer getThickStart() { return thickStart; }
	public Integer getThickEnd() { return thickEnd; }
	public Color getItemRgb() { return itemRgb; }
	public int getBlockCount() { return blockCount; }
	public int[] getBlockSizes() { return blockSizes; }
	public int[] getBlockStarts() { return blockStarts; }
    
    public static Pattern chromPattern = Pattern.compile("chr(.*)");
	
	public BEDLine(String line) { 
		String[] array = line.split("\t");
		chrom = array[0];
        
        Matcher m = chromPattern.matcher(chrom);
        if(m.matches()) { 
            chrom = m.group(1);
        }
        
		chromStart = Integer.parseInt(array[1]);
		chromEnd = Integer.parseInt(array[2]);

		name = array.length > 3 ? array[3] : null;
		score = array.length > 4 ? Double.parseDouble(array[4]) : null;
		strand = array.length > 5 ? array[5].charAt(0) : null;
		thickStart = array.length > 6 ? Integer.parseInt(array[6]) : null;
		thickEnd = array.length > 7 ? Integer.parseInt(array[7]) : null;
		
		String[] aa = array.length > 8 ? array[8].split(",") : null;
		if(array.length > 8) { 
			int r = Integer.parseInt(aa[0]), 
			g = Integer.parseInt(aa[1]), 
			b = Integer.parseInt(aa[2]);
			itemRgb = new Color(r, g, b);
		}
		
		blockCount = array.length > 9 ? Integer.parseInt(array[9]) : null;
		
		aa = array.length > 10 ? array[10].split(",") : null;
		if(array.length > 10) { 
			blockSizes = new int[blockCount];
			for(int i = 0; i < blockCount; i++) { 
				blockSizes[i] = Integer.parseInt(aa[i]);
			}
		}

		aa = array.length > 11 ? array[11].split(",") : null;
		if(array.length > 11) { 
			blockStarts = new int[blockCount];
			for(int i = 0; i < blockCount; i++) { 
				blockStarts[i] = Integer.parseInt(aa[i]);
			}
		}
	}
}

