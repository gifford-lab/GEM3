package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.utils.stats.StatUtil;

public class BEDFileReader extends AlignmentFileReader {

	public BEDFileReader(File f, Genome g, boolean useNonUnique, int idStart) {
		super(f,g,-1,useNonUnique, idStart);
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
	        	if(line.charAt(0)!='#'){
		            String[] words = line.split("\\t");
		            String chr="."; char strand = '.';
		            int start=0, end=0;
		            
		            //String ID = words[3]; //No reliable ID for BED format, so treat EVERY hit as a new/unique read
		            
		            
	            	//if(ID.equals(lastID)){
	            	//	currReadHitCount++;
	            	//}else{
	            		if(currRead!=null){
	            			currRead.setNumHits(currReadHitCount);
	            			//Add the hits to the data structure
	            			addHits(currRead);
	            			currRead=null;
	            		}
	            		currReadHitCount=1;            			
	            	//}
	            	//String tag = words[3];
	            	//if(tag.equals("U") || (useNonUnique && words.length>9 && tag.charAt(0)=='R')){
            			chr = words[0];
            			String[] tmp = chr.split("\\.");
            			chr=tmp[0].replaceFirst("chr", "");
            			chr=chr.replaceFirst("^>", "");
            			start = new Integer(words[1]).intValue();
            			end = new Integer(words[2]).intValue();
            			if(readLength==-1)
    	    				readLength = end-start+1;
            			strand = words[5].charAt(0);
    					ReadHit currHit = new ReadHit(gen,currID,chr, start, end, strand);
    					currID++;
    					//if(!ID.equals(lastID) || currRead==null){
    					currRead = new Read((int)totalWeight);
    					totalWeight++;
    	            	//}
    					currRead.addHit(currHit);
	    			//}
	            	//lastID=ID;
	        	}
            }
	        if(currRead!=null){
    			currRead.setNumHits(currReadHitCount);
    			//Add the hits to the data structure
    			addHits(currRead);
    		}
	        reader.close();
	        populateArrays();
	        
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}//end of countReads method
	
}//end of BEDFileReader class
