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


    /* to do the matching between the two files, we need to scan back and forth between them.
       Since Picard gives us an iterator to read through the input SAM/BAM file, you can think
       of the input stream as the concatenation of, eg, leftbuffer and leftiter.  The current position
       is the position at the head of leftiter behind all the element of leftbuffer.  Storing the previously
       read elements in the arraylist lets us look back through them to try to match them up to 
       a record from the right rights.
    */
    public static ArrayList<SAMRecord> leftbuffer, rightbuffer;
    public static CloseableIterator<SAMRecord> leftiter, rightiter;

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

        leftbuffer = new ArrayList<SAMRecord>();
        rightbuffer = new ArrayList<SAMRecord>();

        leftiter = leftreader.iterator();
        rightiter = rightreader.iterator();

        /* these are the sets of records for the same read that we're
           going to dump to output */
        Collection<SAMRecord> leftrecords = new ArrayList<SAMRecord>();
        Collection<SAMRecord> rightrecords = new ArrayList<SAMRecord>();

        boolean keepgoing = true;
        int added = 0;
        while (keepgoing) {
            SAMRecord left = nextLeft();
            SAMRecord right = nextRight();
            if (left == null || right == null) {
                break;
            }
            leftbuffer.add(left);
            rightbuffer.add(right);            
            /* if the reads from the two files have the same ID, save them
               in the buffer and keep reading.  We do this because the next line
               from one or both files may have the same ID, so we can't output
               anything yet.
            */
            if (left.getReadName().equals(right.getReadName())) {
                continue;
            } 
            /* If there are more reads in both files and we've recently
               run the match-reads-between-files loop, don't run it again yet
            */
            if (leftiter.hasNext() && rightiter.hasNext() && added++ < 1000) {
                continue;
            }           
            added = 0;

            /* this just loops over the left reads trying to match them to right reads.
               If it does find a match at index i in left and index j in right, then
               we can get rid of everything prior to i and j because we know it didn't match
               and we assume the reads are in the same order in both files (even though we
               don't know what that order is*/
            for (int i = 0; i < leftbuffer.size(); i++) {
                int j = 0;
                while (j < rightbuffer.size() && !leftbuffer.get(i).getReadName().equals(rightbuffer.get(j).getReadName())) {
                    j++;
                }
                if (j == rightbuffer.size()) {
                    continue;
                }
                if (debug) {
                    System.err.println(String.format("Found match of %s at %d and %d",leftbuffer.get(i).getReadName(),i,j));
                }

                /* having found a match, find the rest of the reads with that ID  and output */
                int k = i;
                int l = j;
                do {
                    leftrecords.add(leftbuffer.get(k++));
                } while (k < leftbuffer.size() && leftbuffer.get(i).getReadName().equals(leftbuffer.get(k).getReadName()));
                do {
                    rightrecords.add(rightbuffer.get(l++));
                } while (l < rightbuffer.size() && leftbuffer.get(i).getReadName().equals(rightbuffer.get(l).getReadName()));
                dumpRecords(leftrecords, rightrecords);
                leftrecords.clear();
                rightrecords.clear();                
                
                leftbuffer.subList(0,k).clear();
                rightbuffer.subList(0,l).clear();

                i = -1;                
            }
            /* if there's nothing remaining in the left file and we've already tried matching everything in left to
               everything we've read from right, there's no point keeping the right buffer around any more. */
            if (!leftiter.hasNext()) {
                rightbuffer.clear();
            }
            if (!rightiter.hasNext()) {
                leftbuffer.clear();
            }
            if (debug) {
                System.err.println("li.hn " + leftiter.hasNext() + " lb.size " + leftbuffer.size() + 
                                   "ri.hn " + rightiter.hasNext() + " rb.size " + rightbuffer.size());                
            }
            keepgoing = (leftiter.hasNext() || leftbuffer.size() > 0) &&
                (rightiter.hasNext() || rightbuffer.size() > 0);
        }
        dumpRecords(leftrecords, rightrecords);
    }
    /* pull next aligned read from either the input file or, if at EOF, the buffered list */
    public static SAMRecord nextLeft() {
        SAMRecord result = null;
        while (result == null) {
            if (leftiter.hasNext()) {
                result = leftiter.next();
            } else if (leftbuffer.size() > 0) {
                result = leftbuffer.remove(leftbuffer.size() -1);
            } else {
                return null;
            }
            result.setReadName(result.getReadName().replaceAll("/\\d$",""));
            if (result.getReferenceName().equals("*")) {
                result = null;
            }
        }
        return result;
    }
    public static SAMRecord nextRight() {
        SAMRecord result = null;
        while (result == null) {
            if (rightiter.hasNext()) {
                result = rightiter.next();
            } else if (rightbuffer.size() > 0) {
                result = rightbuffer.remove(rightbuffer.size() -1);
            } else {
                return null;
            }
            result.setReadName(result.getReadName().replaceAll("/\\d$",""));
            if (result.getReferenceName().equals("*")) {
                result = null;
            }
        }
        return result;
    }    

    public static void dumpRecords(Collection<SAMRecord> lefts,
                                   Collection<SAMRecord> rights) {
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

    
        
        



}