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

public class ElandFileReader extends AlignmentFileReader {
	public ElandFileReader(File f){
		super(f);
	}
	public ElandFileReader(File f, Genome g, int misMatch, boolean useNonUnique, int idStart,
			HashMap<String, Integer> chrom2ID, HashMap<Integer,String> id2Chrom) {
		super(f,g,misMatch,useNonUnique, idStart, chrom2ID, id2Chrom);
	}

	//Estimate chromosome lengths
	protected void estimateGenome(File f) {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(f));
			String line;
			while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\s+");
		            if(readLength==-1)
	    				readLength = words[1].length();
		            String chr = words[6];
		            String[] tmp = chr.split("\\.");
	            	chr=tmp[0].replaceFirst("chr", "");
	            	chr=chr.replaceFirst("^>", "");
	            	int max = new Integer(words[7]).intValue()+readLength;
	            	
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
	            String[] words = line.split("\\t");
	            String chr="."; char strand = '.';
	            int start=0, end=0;
	            
	            String ID = words[0];
	            if(readLength==-1)
    				readLength = words[1].length();
	            
            	if(ID.equals(lastID)){
            		currReadHitCount++;
            	}else{
            		if(currRead!=null){
            			currRead.setNumHits(currReadHitCount);
            			//Add the hits to the data structure
            			addHits(currRead);
            			currRead=null;
            		}
            		currReadHitCount=1;            			
            	}
            	String tag = words[2];
            	if(tag.charAt(0)=='U' || (useNonUnique && words.length>8 && tag.charAt(0)=='R')){
            		int mis = Integer.parseInt(tag.substring(1, 2));
            		if(mis<=misMatch){
            			chr = words[6];
            			String[] tmp = chr.split("\\.");
            			chr=tmp[0].replaceFirst("^chr", "");
            			start = new Integer(words[7]).intValue();
            			end =start+readLength-1;
            			strand = words[8].equals("F") ? '+' : '-';
    					ReadHit currHit = new ReadHit(gen,currID,chr, start, end, strand, 1, mis);
    					currID++;
    					if(!ID.equals(lastID) || currRead==null){
    						currRead = new Read((int)totalWeight);
    						totalWeight++;
    	            	}currRead.addHit(currHit);
    				}
            	}
            	lastID=ID;
            }
	        if(currRead!=null){
    			currRead.setNumHits(currReadHitCount);
    			//Add the hits to the data structure
    			addHits(currRead);
    		}
	        reader.close();
	        populateArrays();
	        
		} catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of ElandFileReader class
