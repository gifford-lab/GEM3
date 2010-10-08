package edu.mit.csail.cgs.deepseq.utilities;

import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.mit.csail.cgs.utils.RealValuedHistogram;

public class SAMStats {

	private int totalHits=0, LHits=0, RHits=0;
	private int singleEnd=0, properPair=0, unMapped=0, properPairL=0, properPairR=0, pairMapped=0, notPrimary=0; 
	private double weight = 0;
	private int pairedEndSameChr=0, pairedEndDiffChr=0;
	private RealValuedHistogram histo;
	
	public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("i","insdist",false,"print insert size distribution");
        options.addOption("s","stats",false,"print mapping stats");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );       
        
        SAMStats s = new SAMStats();
        if(cl.hasOption("insdist"))
        	s.printInsertDistrib();
        if(cl.hasOption("stats"))
        	s.printStats();
    }
	
	public SAMStats(){
		histo = new RealValuedHistogram(0, 1000, 100);
		SAMFileReader reader = new SAMFileReader(System.in);
		reader.setValidationStringency(ValidationStringency.SILENT);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            processRecord(record);
        }
        iter.close();
        reader.close();        
    }
	
	public void processRecord(SAMRecord r){
		if(r.getMateUnmappedFlag())
			unMapped++;
		else{
			totalHits++;
			weight += 1/(float)r.getIntegerAttribute("NH");
			
			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag())
					singleEnd++;
				else if(r.getMateReferenceName().equals(r.getReferenceName()))
					pairedEndSameChr++;
				else
					pairedEndDiffChr++;
			}else{
				singleEnd++;
			}
			
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					if(r.getProperPairFlag()){
						properPairL++;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
					if(r.getProperPairFlag())
						properPairR++;
					if(!r.getMateUnmappedFlag()){
						pairMapped++;
					}
				}
			}else{
				notPrimary++;
			}
		}
	}
	
	
	public void printInsertDistrib(){
		histo.printContents();
	}

	public void printStats(){
		System.out.println("\nTotal Mappings:\t"+totalHits);
		System.out.println("Total Weight:\t"+weight);
		System.out.println("Left Mappings\t"+LHits);
		System.out.println("Right Mappings\t"+RHits);
		System.out.println("Paired Hit Mappings\t"+pairMapped);
		System.out.println("Proper pairs (L):\t"+properPairL);
		System.out.println("Proper pairs (R):\t"+properPairR);
		System.out.println("UnMapped:\t"+unMapped);
		System.out.println("NotPrimary:\t"+notPrimary);
		System.out.println("SingleEnd Mappings:\t"+singleEnd);
		System.out.println("PairedEnd Mappings (same chr):\t"+pairedEndSameChr);
		System.out.println("PairedEnd Mappings (diff chr):\t"+pairedEndDiffChr);
	}
}
