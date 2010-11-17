package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

/**
 * Reads SAM or BAM data on stdin.
 * Produces a file on stdout in the format expected by ImportHits.
 * The weight for a hit is 1/(# of hits for that read)
 * 
 * Options:	--nosuboptimal (flag to only take the hits with the minimum number of mismatches)
 * 			--uniquehits (flag to only print 1:1 read to hit mappings)
 * 
 * nosuboptimal is applied before uniquehits
 */

public class SAMToReadDB {

    public static boolean uniqueOnly;
    public static boolean filterSubOpt;

    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("s","nosuboptimal",false,"do not include hits whose score is not equal to the best score for the read");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	uniqueOnly = cl.hasOption("uniquehits");
    	filterSubOpt = cl.hasOption("nosuboptimal");
        String line;
        String lastRead = "";        
        SAMFileReader reader = new SAMFileReader(System.in);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        Collection<SAMRecord> byRead = new ArrayList<SAMRecord>();
        String lastread = null;
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if (record.getReadUnmappedFlag()) {continue; }
            if (lastread == null || !lastread.equals(record.getReadName())) {
                dumpRecords(byRead);
                byRead.clear();
            }
            lastread = record.getReadName();
            byRead.add(record);
            
        }
        dumpRecords(byRead);
        iter.close();
        reader.close();
    }       
    public static Collection<SAMRecord> filterNoChrom(Collection<SAMRecord> input) {
        if (input.size() == 0) {
            return input;
        }
        Collection<SAMRecord> output = new ArrayList<SAMRecord>();
        for (SAMRecord r : input) {
            if (!r.getReferenceName().equals("*")) {
                output.add(r);
            }
        }
        return output;
    }

    public static Collection<SAMRecord> filterSubOpt(Collection<SAMRecord> input) {
        if (input == null || input.size() < 2) {
            return input;
        }
        int maxqual = SAMRecord.NO_MAPPING_QUALITY;
        for (SAMRecord r : input) {
            if (r.getMappingQuality() > maxqual) {
                maxqual = r.getMappingQuality();
            }
        }
        Collection<SAMRecord> output = new ArrayList<SAMRecord>();
        for (SAMRecord r : input) {
            if (maxqual == r.getMappingQuality()) {
                output.add(r);
            }
        }
        return output;
    }
    public static void dumpRecords(Collection<SAMRecord> records) {
        
        int mapcount = records.size();
        if (mapcount == 0) {
            return;
        }
        if (filterSubOpt) {
            records = filterSubOpt(records);
        }        
        if (uniqueOnly && mapcount > 1) {
            return;
        }

        float weight = 1 / ((float)mapcount);

        for (SAMRecord record : records) {
            System.out.println(String.format("%s\t%d\t%s\t%d\t%f",
                                             record.getReferenceName(),
                                             record.getReadNegativeStrandFlag() ? 
                                             record.getAlignmentEnd() : 
                                             record.getAlignmentStart(),
                                             record.getReadNegativeStrandFlag() ? "-" : "+",
                                             record.getReadLength(),
                                             weight));
        }                                             
    }
}