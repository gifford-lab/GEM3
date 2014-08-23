package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

public class SAMReader extends AlignmentFileReader{
	public SAMReader(File f){
		super(f);
	}

    public SAMReader(File f, Genome g, int mis, boolean nonUnique, int idSeed,
			HashMap<String, Integer> chrom2ID, HashMap<Integer,String> id2Chrom) {
	super(f, g, mis, nonUnique, idSeed, chrom2ID, id2Chrom);
    }
    
	protected void estimateGenome(File f) {
		HashMap<String, Integer> chrLenMap = new HashMap<String, Integer>();
		SAMFileReader reader = new SAMFileReader(f);
		reader.setValidationStringency(ValidationStringency.SILENT);
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
		reader.close();
	}
	
    //Return the total reads and weight
    protected void countReads() {
		readLength=-1;
		totalHits=0;
		totalWeight=0;
		
		SAMFileReader reader = new SAMFileReader(inFile);
		reader.setValidationStringency(ValidationStringency.LENIENT);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
		String lastread = null;
		while (iter.hasNext()) {
		    currID++;
		    SAMRecord record = iter.next();
		    if(readLength ==-1)
		    	readLength = record.getReadLength();
		    
		    if (record.getReadUnmappedFlag()) {
		    	continue; 
		    }
		    if (lastread == null || !lastread.equals(record.getReadName())) {
		    	processRead(byRead);
		    	byRead.clear();
		    }
		    lastread = record.getReadName();
		    byRead.add(record);
			    
		}
		processRead(byRead);
		iter.close();
		reader.close();
		populateArrays();
    }//end of countReads method
    
    protected void processRead(Collection<SAMRecord> records) {
        int mapcount = records.size();
        if (mapcount == 0) {
            return;
        }
//        if (!useNonUnique && mapcount > 1) {		// old code: this will ignore paired-end reads
//            return;
//        }
        
        boolean paired = false;
    	for (SAMRecord record : records) {
    		if (record.getReadPairedFlag()){
    			paired = true;
    			break;
    		}
    	}        
    	
    	if (!useNonUnique) {		// accepting only uniquely-mapped reads 
        	if (!paired && mapcount > 1)
        		return;
        }
        
        float weight = 1 / ((float)mapcount);
        if (paired)
        	weight = 1;
        Read currRead = new Read((int)totalWeight);
		for (SAMRecord record : records) {
		    int start = record.getAlignmentStart();
		    int end = record.getAlignmentEnd();
		    ReadHit currHit = new ReadHit(gen,
						  currID,
						  record.getReferenceName().replaceFirst("^chr", ""), 
						  start, end, 
						  record.getReadNegativeStrandFlag() ? '-' : '+',
						  weight);
		    currRead.addHit(currHit);
		    currID++;
		}
		currRead.setNumHits(mapcount);
		addHits(currRead);
		if (paired)
			totalWeight+=mapcount;
		else
			totalWeight++;
    }//end of processRead

}
