package edu.mit.csail.cgs.projects.readdb;

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

	private int totalHits=0, LHits=0, RHits=0, totalHitBP=0, LHitBP=0, RHitBP=0;
	private int singleEnd=0, properPair=0, unMapped=0, properPairL=0, properPairR=0, pairMapped=0, notPrimary=0;
	private int singleEndBP=0, properPairBP=0, unMappedBP=0, properPairLBP=0, properPairRBP=0, pairMappedBP=0, notPrimaryBP=0;
	private int uniquelyMapped=0, uniquelyMappedBP=0;
	private double weight=0, weightBP=0;;
	private int pairedEndSameChr=0, pairedEndDiffChr=0;
	private int pairedEndSameChrBP=0, pairedEndDiffChrBP=0;
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
		if(r.getReadUnmappedFlag())
			unMapped++;
		else{
			totalHits++;
			int len = r.getReadLength();
			totalHitBP+=len;
			int count = r.getIntegerAttribute("NH");
			if(count==1){
				uniquelyMapped++;
				uniquelyMappedBP+=len;
			}
			
			weight += 1/(float)count;
			weightBP += (1/(float)count)*len;
			
			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag()){
					singleEnd++;
					singleEndBP+=len;
				}else{
					pairMapped++;
					pairMappedBP+=len;
					if(r.getMateReferenceName().equals(r.getReferenceName())){
						pairedEndSameChr++;
						pairedEndSameChrBP+=len;
					}else{
						pairedEndDiffChr++;
						pairedEndDiffChrBP+=len;
					}
				}
			}else{
				singleEnd++;
				singleEndBP+=len;
			}
			
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					LHitBP+=len;
					if(r.getReadPairedFlag() && r.getProperPairFlag()){
						properPairL++;
						properPairLBP+=len;
						properPair++;
						properPairBP+=len;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
					RHitBP+=len;
					if(r.getProperPairFlag()){
						properPairR++;
						properPairRBP+=len;
						properPair++;
						properPairBP+=len;
					}
				}
			}else{
				notPrimary++;
				notPrimaryBP+=len;
			}
		}
	}
	
	
	public void printInsertDistrib(){
		histo.printContents();
	}

	public void printStats(){
		System.out.println("\nTotalHits:\t"+totalHits+"\t"+totalHitBP+" bp");
		System.out.println("LeftHits\t"+LHits+"\t"+LHitBP+" bp");
		System.out.println("RightHits\t"+RHits+"\t"+RHitBP+" bp");
		System.out.println("MappedSeq:\t"+(int)weight+"\t"+weightBP+" bp");
		System.out.println("UniquelyMapped:\t"+uniquelyMapped+"\t"+uniquelyMappedBP+" bp");
		int nonU = (int)weight - uniquelyMapped;
		int nonUBP = (int)weightBP - uniquelyMappedBP;
		System.out.println("NonUniquelyMapped:\t"+nonU+"\t"+nonUBP+" bp");
		System.out.println("SingleEndMapped:\t"+singleEnd+"\t"+singleEndBP+" bp");
		System.out.println("PairedEndMapped\t"+pairMapped+"\t"+pairMappedBP+" bp");
		System.out.println("ProperPairs:\t"+properPair+"\t"+properPairBP+" bp");
		System.out.println("ProperPairsL:\t"+properPairL+"\t"+properPairLBP+" bp");
		System.out.println("ProperPairsR:\t"+properPairR+"\t"+properPairRBP+" bp");
		System.out.println("PairedEndMapped_SameChr:\t"+pairedEndSameChr+"\t"+pairedEndSameChrBP+" bp");
		System.out.println("PairedEndMapped_DiffChr:\t"+pairedEndDiffChr+"\t"+pairedEndDiffChrBP+" bp");
		System.out.println("UnMapped:\t"+unMapped+"\t"+unMappedBP+" bp");
		System.out.println("NotPrimary:\t"+notPrimary+"\t"+notPrimaryBP+" bp");
	}
}
