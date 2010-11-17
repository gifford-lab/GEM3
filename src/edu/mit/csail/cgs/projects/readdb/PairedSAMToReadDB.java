package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;


/**
 * Reads two files of SAM or BAM data and produces output on stdout in the
 * format expected by ImportHits.  Both files must be sorted in the same order.
 * Only reads present in both files will be included in the output (on stdout).
 * 
 * The matching of reads between files is done by stripping "/\d" from the end of the 
 * read name, as reads usually end in /1 or /2.
 *
 * Usage:
 * java PairedSAMToReadDB --left leftreads.bam --right rightreads.bam
 *
 *
 * Options:	--nosuboptimal (flag to only take the hits with the minimum number of mismatches)
 * 			--uniquehits (flag to only print 1:1 read to hit mappings)
 * 
 * nosuboptimal is applied before uniquehits
 *
 * Output columns are
 * 1) left chromname
 * 2) left position
 * 3) left strand
 * 4) left readlen
 * 5) right chromname
 * 6) right position
 * 7) right strand
 * 8) right length
 * 9) weight
 */


public class PairedSAMToReadDB {

    public static boolean uniqueOnly, filterSubOpt, debug;
    public static int chunksize = 2000;
    public static ArrayList<SAMRecord> leftbuffer, rightbuffer;
    public static CloseableIterator<SAMRecord> leftiter, rightiter;

    public static void dumpRecords(Collection<SAMRecord> lefts,
                                   Collection<SAMRecord> rights) {
        if (filterSubOpt) {
            lefts = SAMToReadDB.filterSubOpt(SAMToReadDB.filterNoChrom(lefts));
            rights = SAMToReadDB.filterSubOpt(SAMToReadDB.filterNoChrom(rights));
        }
        int mapcount = lefts.size() * rights.size();
        if (mapcount == 0) {
            return;
        }
        if (uniqueOnly && mapcount > 1) {
            return;
        }
        float weight = 1 / ((float)mapcount);
        for (SAMRecord left : lefts) {
            for (SAMRecord right : rights) {
                System.out.println(String.format("%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%f",
                                                 left.getReferenceName(),
                                                 left.getReadNegativeStrandFlag() ? 
                                                 left.getAlignmentEnd() : 
                                                 left.getAlignmentStart(),
                                                 left.getReadNegativeStrandFlag() ? "-" : "+",
                                                 left.getReadLength(),
                                                 
                                                 right.getReferenceName(),
                                                 right.getReadNegativeStrandFlag() ? 
                                                 right.getAlignmentEnd() : 
                                                 right.getAlignmentStart(),
                                                 right.getReadNegativeStrandFlag() ? "-" : "+",
                                                 right.getReadLength(),

                                                 weight));                
            }
        }
    }

    public static boolean fillLeft(int n) {
        boolean filled = false;
        for (int i = 0; i < n; i++) {
            if (leftiter.hasNext()) {
                SAMRecord r = leftiter.next();
                r.setReadName(r.getReadName().replaceAll("/\\d$",""));
                leftbuffer.add(r);
                filled = true;
            }
        }
        return filled;
    }
    public static boolean fillRight(int n) {
        boolean filled = false;
        for (int i = 0; i < n; i++) {
            if (rightiter.hasNext()) {
                SAMRecord r = rightiter.next();
                r.setReadName(r.getReadName().replaceAll("/\\d$",""));
                rightbuffer.add(r);
                filled = true;
            }
        }
        return filled;
    }

    public static boolean fill(int n) {
        return fillLeft(n) && fillRight(n);
    }

    public static void makePairs() {
        ArrayList<SAMRecord> leftrecords = new ArrayList<SAMRecord>();
        ArrayList<SAMRecord> rightrecords = new ArrayList<SAMRecord>();
        int loop = 0;
        while (fill(chunksize) || leftbuffer.size() > 0) {
            int clearL = 0, clearR = 0;
            int lbs = Math.min(leftbuffer.size(),10000);
            for (int i = 0; i < lbs; i++) {
                String readname = leftbuffer.get(i).getReadName();
                String nextreadname = i < leftbuffer.size() - 1 ? leftbuffer.get(i+1).getReadName() : null;
                int j = clearR;                
                while (j < rightbuffer.size()) {
                    if (readname.equals(rightbuffer.get(j).getReadName())) {
                        /* having found a match, find the rest of the reads with that ID  and output */
                        int k = i;
                        int l = j;
                        /* make sure we'll be able to find all reads with this name */
                        while (readname.equals(leftbuffer.get(leftbuffer.size()-1).getReadName()) &&
                               fillLeft(chunksize)) {
                        }
                        while (readname.equals(rightbuffer.get(rightbuffer.size()-1).getReadName()) &&
                               fillRight(chunksize)) {
                        }
                        do {
                            leftrecords.add(leftbuffer.get(k++));
                        } while (k < leftbuffer.size() && readname.equals(leftbuffer.get(k).getReadName()));
                        do {
                            rightrecords.add(rightbuffer.get(l++));
                        } while (l < rightbuffer.size() && readname.equals(rightbuffer.get(l).getReadName()));
                        dumpRecords(leftrecords, rightrecords);
                        leftrecords.clear();
                        rightrecords.clear();                
                        clearL = k;
                        clearR = l;
                        i = k-1;                            
                        break;
                    } else if (nextreadname != null && nextreadname.equals(rightbuffer.get(j).getReadName())) {
                        clearR = j;
                        break;
                    }
                    j++;
                }
            }
            if (rightiter.hasNext()) {
                leftbuffer.subList(0,clearL).clear();
                rightbuffer.subList(0,clearR).clear();
            } else {
                leftbuffer.clear();
            }
            System.err.println(String.format("loop %d, left %d, right %d",
                                             loop++, leftbuffer.size(), rightbuffer.size()));
        }
    }

    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("l","left",true,"filename of left side of read");
        options.addOption("r","right",true,"filename of right side of read");
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("s","nosuboptimal",false,"do not include hits whose score is not equal to the best score for the read");
        options.addOption("D","debug",false,"enable debugging spew?");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	uniqueOnly = cl.hasOption("uniquehits");
    	filterSubOpt = cl.hasOption("nosuboptimal");
        debug = cl.hasOption("debug");
        String leftfile = cl.getOptionValue("left");
        String rightfile = cl.getOptionValue("right");

        SAMFileReader leftreader = new SAMFileReader(new FileInputStream(leftfile));
        SAMFileReader rightreader = new SAMFileReader(new FileInputStream(rightfile));
        leftiter = leftreader.iterator();
        rightiter = rightreader.iterator();

        leftbuffer = new ArrayList<SAMRecord>();
        rightbuffer = new ArrayList<SAMRecord>();

        makePairs();
    }
}