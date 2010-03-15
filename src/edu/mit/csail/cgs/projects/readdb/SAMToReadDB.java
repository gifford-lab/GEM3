package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;
import net.sf.samtools.*;

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

    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("s","nosuboptimal",false,"do not include hits whose score is not equal to the best score for the read");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	boolean uniqueOnly = cl.hasOption("uniquehits");
    	boolean filterSubOpt = cl.hasOption("nosuboptimal");
        String line;
        String lastRead = "";        
        SAMFileReader reader = new SAMFileReader(System.in);
        Iterator<SAMRecord> iter = reader.iterator();
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if (record.getReadUnmappedFlag()) {continue; }
            int mapcount = (Integer)record.getAttribute("NH");
            
            if (uniqueOnly && mapcount > 1) {
                continue;
            }
            float weight = 1 / ((float)mapcount);

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