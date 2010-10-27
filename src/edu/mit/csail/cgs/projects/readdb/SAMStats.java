package edu.mit.csail.cgs.projects.readdb;

import java.io.IOException;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
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
	private double singleEndBP=0, properPairBP=0, unMappedBP=0, properPairLBP=0, properPairRBP=0, pairMappedBP=0, notPrimaryBP=0;
	private int uniquelyMapped=0, uniquelyMappedBP=0;
	private int junctions=0, junctionsBP=0;
	private double weight=0, weightBP=0;;
	private int pairedEndSameChr=0, pairedEndDiffChr=0;
	private double pairedEndSameChrBP=0, pairedEndDiffChrBP=0;
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
			double dlen = (double)len;
			totalHitBP+=dlen;
			int count = r.getIntegerAttribute("NH");
			if(count==1){
				uniquelyMapped++;
				uniquelyMappedBP+=dlen;
			}
			
			weight += 1/(float)count;
			weightBP += (1/(float)count)*dlen;
			
			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag()){
					singleEnd++;
					singleEndBP+=dlen;
				}else{
					pairMapped++;
					pairMappedBP+=dlen;
					if(r.getMateReferenceName().equals(r.getReferenceName())){
						pairedEndSameChr++;
						pairedEndSameChrBP+=dlen;
					}else{
						pairedEndDiffChr++;
						pairedEndDiffChrBP+=dlen;
					}
				}
			}else{
				singleEnd++;
				singleEndBP+=dlen;
			}
			
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
    		if(blocks.size()==2){
    			junctions++;
    			junctionsBP+=dlen;
    		}
			
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					LHitBP+=len;
					if(r.getReadPairedFlag() && r.getProperPairFlag()){
						properPairL++;
						properPairLBP+=dlen;
						properPair++;
						properPairBP+=dlen;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
					RHitBP+=dlen;
					if(r.getProperPairFlag()){
						properPairR++;
						properPairRBP+=dlen;
						properPair++;
						properPairBP+=dlen;
					}
				}
			}else{
				notPrimary++;
				notPrimaryBP+=dlen;
			}
		}
	}
	
	
	public void printInsertDistrib(){
		histo.printContents();
	}

	public void printStats(){
		System.out.println("\nTotalHits:\t"+totalHits+"\t"+String.format("%.0f",totalHitBP)+" bp");
		System.out.println("LeftHits:\t"+LHits+"\t"+String.format("%.0f",LHitBP)+" bp");
		System.out.println("RightHits:\t"+RHits+"\t"+String.format("%.0f",RHitBP)+" bp");
		System.out.println("MappedSeq:\t"+(int)weight+"\t"+weightBP+" bp");
		System.out.println("UniquelyMapped:\t"+uniquelyMapped+"\t"+String.format("%.0f",uniquelyMappedBP)+" bp");
		int nonU = (int)weight - uniquelyMapped;
		int nonUBP = (int)weightBP - uniquelyMappedBP;
		System.out.println("NonUniquelyMapped:\t"+nonU+"\t"+String.format("%.0f",nonUBP)+" bp");
		System.out.println("SingleEndMapped:\t"+singleEnd+"\t"+String.format("%.0f",singleEndBP)+" bp");
		System.out.println("Junctions:\t"+junctions+"\t"+String.format("%.0f",junctionsBP)+" bp");
		System.out.println("PairedEndMapped:\t"+pairMapped+"\t"+String.format("%.0f",pairMappedBP)+" bp");
		System.out.println("ProperPairs:\t"+properPair+"\t"+String.format("%.0f",properPairBP)+" bp");
		System.out.println("ProperPairsL:\t"+properPairL+"\t"+String.format("%.0f",properPairLBP)+" bp");
		System.out.println("ProperPairsR:\t"+properPairR+"\t"+String.format("%.0f",properPairRBP)+" bp");
		System.out.println("PairedEndMapped_SameChr:\t"+pairedEndSameChr+"\t"+String.format("%.0f",pairedEndSameChrBP)+" bp");
		System.out.println("PairedEndMapped_DiffChr:\t"+pairedEndDiffChr+"\t"+String.format("%.0f",pairedEndDiffChrBP)+" bp");
		System.out.println("UnMapped:\t"+unMapped+"\t"+String.format("%.0f",unMappedBP)+" bp");
		System.out.println("NotPrimary:\t"+notPrimary+"\t"+String.format("%.0f",notPrimaryBP)+" bp");
	}
}
