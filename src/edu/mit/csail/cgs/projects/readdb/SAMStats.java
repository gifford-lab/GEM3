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

	private Double totalHits=0.0, LHits=0.0, RHits=0.0, totalHitBP=0.0, LHitBP=0.0, RHitBP=0.0;
	private Double singleEnd=0.0, properPair=0.0, unMapped=0.0, properPairL=0.0, properPairR=0.0, pairMapped=0.0, notPrimary=0.0;
	private Double singleEndBP=0.0, properPairBP=0.0, unMappedBP=0.0, properPairLBP=0.0, properPairRBP=0.0, pairMappedBP=0.0, notPrimaryBP=0.0;
	private Double uniquelyMapped=0.0, uniquelyMappedBP=0.0;
	private Double junctions=0.0, junctionsBP=0.0;
	private Double weight=0.0, weightBP=0.0;;
	private Double pairedEndSameChr=0.0, pairedEndDiffChr=0.0;
	private Double pairedEndSameChrBP=0.0, pairedEndDiffChrBP=0.0;
	private Boolean bowtie2=false;
	private RealValuedHistogram histo;
	
	public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("i","insdist",false,"print insert size distribution");
        options.addOption("s","stats",false,"print mapping stats");
        options.addOption("bt2",false,"input is from bowtie2");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );       
        
        SAMStats s = new SAMStats(cl.hasOption("bt2"));
        if(cl.hasOption("insdist"))
        	s.printInsertDistrib();
        if(cl.hasOption("stats"))
        	s.printStats();
    }
	
	public SAMStats(boolean bowtie2){
		this.bowtie2 = bowtie2;
		histo = new RealValuedHistogram(0, 1000, 100);
		SAMFileReader reader = new SAMFileReader(System.in);
		reader.setValidationStringency(ValidationStringency.SILENT);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if(bowtie2)
            	processBT2SAMRecord(record);
            else
            	processSAMRecord(record);
        }
        iter.close();
        reader.close();        
    }
	
	public void processSAMRecord(SAMRecord r){
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
	
	public void processBT2SAMRecord(SAMRecord r){
		if(r.getReadUnmappedFlag())
			unMapped++;
		else{
			totalHits++;
			int len = r.getReadLength();
			double dlen = (double)len;
			totalHitBP+=dlen;
			int count =1; //TODO: Fix this if using bowtie2 for multi-mapping reads
			
			int primAScore = r.getIntegerAttribute("AS");
	    	int secAScore=-1000000;
	    	if(r.getIntegerAttribute("XS")!=null)
	    		secAScore = r.getIntegerAttribute("XS");
	    	boolean currUnique = primAScore > secAScore ? true : false;
	    	
			if(count==1 && currUnique){
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
		System.out.println("\nTotalHits:\t"+String.format("%.0f",+totalHits)+"\t"+String.format("%.0f",totalHitBP)+" bp");
		System.out.println("LeftHits:\t"+String.format("%.0f",LHits)+"\t"+String.format("%.0f",LHitBP)+" bp");
		System.out.println("RightHits:\t"+String.format("%.0f",RHits)+"\t"+String.format("%.0f",RHitBP)+" bp");
		System.out.println("MappedSeq:\t"+String.format("%.0f",weight)+"\t"+String.format("%.0f",weightBP)+" bp");
		System.out.println("UniquelyMapped:\t"+String.format("%.0f",uniquelyMapped)+"\t"+String.format("%.0f",uniquelyMappedBP)+" bp");
		double nonU = weight - uniquelyMapped;
		double nonUBP = weightBP - uniquelyMappedBP;
		System.out.println("NonUniquelyMapped:\t"+String.format("%.0f",nonU)+"\t"+String.format("%.0f",nonUBP)+" bp");
		System.out.println("SingleEndMapped:\t"+String.format("%.0f",singleEnd)+"\t"+String.format("%.0f",singleEndBP)+" bp");
		System.out.println("PairedEndMapped:\t"+String.format("%.0f",pairMapped)+"\t"+String.format("%.0f",pairMappedBP)+" bp");
		System.out.println("ProperPairs:\t"+String.format("%.0f",properPair)+"\t"+String.format("%.0f",properPairBP)+" bp");
		System.out.println("ProperPairsL:\t"+String.format("%.0f",properPairL)+"\t"+String.format("%.0f",properPairLBP)+" bp");
		System.out.println("ProperPairsR:\t"+String.format("%.0f",properPairR)+"\t"+String.format("%.0f",properPairRBP)+" bp");
		System.out.println("PairedEndMapped_SameChr:\t"+String.format("%.0f",pairedEndSameChr)+"\t"+String.format("%.0f",pairedEndSameChrBP)+" bp");
		System.out.println("PairedEndMapped_DiffChr:\t"+String.format("%.0f",pairedEndDiffChr)+"\t"+String.format("%.0f",pairedEndDiffChrBP)+" bp");
		System.out.println("UnMapped:\t"+String.format("%.0f",unMapped)+"\t"+String.format("%.0f",unMappedBP)+" bp");
		System.out.println("NotPrimary:\t"+String.format("%.0f",notPrimary)+"\t"+String.format("%.0f",notPrimaryBP)+" bp");
		if(!bowtie2){
			System.out.println("Junctions:\t"+String.format("%.0f",junctions)+"\t"+String.format("%.0f",junctionsBP)+" bp");
		}
	}
}
