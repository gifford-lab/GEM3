package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

/**
 * Here:  <br>
 * <tt>totalHits</tt>: total number of hits for this file. (Only unique reads are
 * counted if the corresponding boolean parameter <tt>useNonUnique</tt> is set to false. <br>
 * <tt>totalWeight = totalHits</tt>
 * @author shaunmahony
 *
 */
public class BowtieFileReader extends AlignmentFileReader {

	public BowtieFileReader(File f, Genome g, boolean useNonUnique, int idStart,
			HashMap<String, Integer> chrom2ID, HashMap<Integer,String> id2Chrom) {
		super(f,g,-1,useNonUnique, idStart, chrom2ID, id2Chrom);
	}	
		
	//Estimate chromosome lengths
	protected void estimateGenome() {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			String line;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\s+");
		            if(readLength==-1)
	    				readLength = words[4].length();
		            String chr = words[2];
		            String[] tmp = chr.split("\\.");
	            	chr=tmp[0].replaceFirst("chr", "");
	            	chr=chr.replaceFirst("^>", "");
	            	int max = new Integer(words[3]).intValue()+readLength;
	            	
	        		if(!chrLenMap.containsKey(chr) || chrLenMap.get(chr)<max)
						chrLenMap.put(chr, max);
	        	}
			}
			gen=new Genome("Genome", chrLenMap);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	//Return the total reads and weight
	protected void countReads() {
		try {
			readLength=-1;
			totalHits=0;
			totalWeight=0;
			BufferedReader reader = new BufferedReader(new FileReader(inFile));
			String line, lastID="";
			double currReadHitCount=0;
			Read currRead=null;
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if (line.length()==0) continue;
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr="."; char strand = '.';
		            int start=0, end=0;
		            
		            String ID = words[0];
		            if(readLength==-1)
	    				readLength = words[4].length();
		            
		            boolean newRead=true;
	            	if(ID.equals(lastID)){
	            		currReadHitCount++;
	            		newRead=false;
	            	}else{
	            		if(currRead!=null){
	            			Read filterRead = currReadHitCount==1 ? currRead : currRead.filter();
	            			if(filterRead.getNumHits()==1 || useNonUnique){
	        	        		//Add the hits to the data structure
	        	        		addHits(filterRead);
	        	        		totalWeight++;
	        	        	}currRead=null;
	            		}
	            		currReadHitCount=1;            			
	            	}
	            	int mis=0;
	            	if(words.length>7 && words[7].length()>1){
	            		mis = words[7].split(",").length;
	            	}
	            	chr = words[2];
        			String[] tmp = chr.split("\\.");
        			chr=tmp[0].replaceFirst("chr", "");
        			chr=chr.replaceFirst("^>", "");
        			start = new Integer(words[3]).intValue();
        			end =start+readLength-1;
        			strand = words[1].charAt(0);
					ReadHit currHit = new ReadHit(gen,currID,chr, start, end, strand, 1, mis);
					currID++;
					if(newRead || currRead==null){
						currRead = new Read((int)totalWeight);						
	            	}
					currRead.addHit(currHit);
    			
	            	lastID=ID;
	        	}
            }
	        if(currRead!=null){
	        	Read filterRead = currReadHitCount==1 ? currRead : currRead.filter();
    			if(filterRead.getNumHits()==1 || useNonUnique){
	        		//Add the hits to the data structure
	        		addHits(filterRead);
	        		totalWeight++;
	        	}
    		}
	        reader.close();
	        populateArrays();
	 
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of BowtieFileReader class
