package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.Read;
import edu.mit.csail.cgs.deepseq.ReadHit;

public class SAMReader extends AlignmentFileReader{

	protected List<SAMReader> readers=new ArrayList<SAMReader>();
	
	public SAMReader(File f, Genome g, int mis, boolean nonUnique, int idSeed) {
		super(f, g, mis, nonUnique, idSeed);
		
	}
	
	//Return the total reads and weight
	protected void countReads() {
		readLength=-1;
		totalHits=0;
		totalWeight=0;
	
		SAMFileReader reader = new SAMFileReader(System.in);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
        String lastread = null;
        while (iter.hasNext()) {
        	currID++;
            SAMRecord record = iter.next();
            if(readLength ==-1)
            	readLength = record.getReadLength();
            
            if (record.getReadUnmappedFlag()) {continue; }
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
        if (!useNonUnique && mapcount > 1) {
            return;
        }

        float weight = 1 / ((float)mapcount);
        Read currRead = new Read((int)totalWeight);
        for (SAMRecord record : records) {
			int start = record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart();
			int end = start + record.getReadLength()-1;
			ReadHit currHit = new ReadHit(gen,
					currID,
					record.getReferenceName(), 
					start, end, 
					record.getReadNegativeStrandFlag() ? '-' : '+');
			currID++;
		}currRead.setNumHits(mapcount);
		totalWeight++;
    }//end of processRead
}
