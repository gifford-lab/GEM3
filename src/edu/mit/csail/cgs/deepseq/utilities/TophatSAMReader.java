package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

public class TophatSAMReader extends AlignmentFileReader{

    public TophatSAMReader(File f, Genome g, int mis, boolean nonUnique, int idSeed) {
    	super(f, g, mis, nonUnique, idSeed);
    }
    
	protected void estimateGenome() {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		SAMFileReader reader = new SAMFileReader(inFile);
		SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
		if(dictionary !=null){
			for(SAMSequenceRecord record : dictionary.getSequences()){
			    String chr = record.getSequenceName().replaceFirst("^chr", "");
			    chrLenMap.put(chr, record.getSequenceLength());
			}
		}else{
			CloseableIterator<SAMRecord> iter = reader.iterator();
			while (iter.hasNext()) {
			    currID++;
			    SAMRecord record = iter.next();
			    String chr = record.getReferenceName().replaceFirst("^chr", "");
			    int max = Math.max(record.getAlignmentEnd(), record.getAlignmentStart());
			    if(!chrLenMap.containsKey(chr) || chrLenMap.get(chr)<max)
					chrLenMap.put(chr, max);
			}
		}
		gen=new Genome("Genome", chrLenMap);
	}
	
    //Return the total reads and weight
    protected void countReads() {
		readLength=-1;
		totalHits=0;
		totalWeight=0;
		
		SAMFileReader reader = new SAMFileReader(inFile);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
		    currID++;
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    float weight = 1/record.getFloatAttribute("NH");
		    if(readLength ==-1)
		    	readLength = record.getReadLength();
		    
	        Read currRead = new Read((int)totalWeight);
	        
	        List<AlignmentBlock> blocks = record.getAlignmentBlocks();
		    for(int a=0; a<blocks.size(); a++){ //Iterate over alignment blocks
		    	AlignmentBlock currBlock = blocks.get(a);
		    	int aStart = currBlock.getReferenceStart();
		    	int aEnd = aStart + currBlock.getLength()-1;
		    	int aLen = currBlock.getLength();
		    	boolean nearbyBlocks=true;
		    	while(nearbyBlocks && a<blocks.size()){
		    		if(blocks.get(a+1).getReferenceStart() - currBlock.getReferenceStart() < record.getReadLength()){
		    			aEnd = blocks.get(a+1).getReferenceStart() + blocks.get(a+1).getLength()-1;
		    			aLen += blocks.get(a+1).getLength();
		    			a++;
		    		}else{
		    			nearbyBlocks=false;
		    		}
		    	}
		    	
		    	ReadHit currHit = new ReadHit(gen,
					  currID,
					  record.getReferenceName().replaceFirst("^chr", ""), 
					  aStart, aEnd, 
					  record.getReadNegativeStrandFlag() ? '-' : '+',
					  weight);
		   
		    	currRead.addHit(currHit);
		    	currID++;
			}	
		    addHits(currRead);
		    totalWeight++;
		}
		iter.close();
		reader.close();
		populateArrays();
    }//end of countReads method
    
}
