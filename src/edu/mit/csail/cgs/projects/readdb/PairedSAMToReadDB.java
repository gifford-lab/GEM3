package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;


/**
 * Reads two files of SAM or BAM data and produces output on stdout in the
 * format expected by ImportHits.  Both files must be sorted in the same order.
 * Only reads present in both files will be included in the output.
 *
 * Options:	--nosuboptimal (flag to only take the hits with the minimum number of mismatches)
 * 			--uniquehits (flag to only print 1:1 read to hit mappings)
 * 
 * nosuboptimal is applied before uniquehits
 */


public class PairedSAMToReadDB {

    public static boolean uniqueOnly;
    public static boolean filterSubOpt;
    public static ArrayList<SAMRecord> leftbuffer, rightbuffer;
    public static CloseableIterator<SAMRecord> leftiter, rightiter;

    public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("l","left",true,"filename of left side of read");
        options.addOption("r","right",true,"filename of right side of read");
        options.addOption("u","uniquehits",false,"only output hits with a single mapping");
        options.addOption("s","nosuboptimal",false,"do not include hits whose score is not equal to the best score for the read");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );            
    	uniqueOnly = cl.hasOption("uniquehits");
    	filterSubOpt = cl.hasOption("nosuboptimal");
        String leftfile = cl.getOptionValue("left");
        String rightfile = cl.getOptionValue("right");

        SAMFileReader leftreader = new SAMFileReader(new FileInputStream(leftfile));
        SAMFileReader rightreader = new SAMFileReader(new FileInputStream(rightfile));

        leftbuffer = new ArrayList<SAMRecord>();
        rightbuffer = new ArrayList<SAMRecord>();

        leftiter = leftreader.iterator();
        rightiter = rightreader.iterator();

        Collection<SAMRecord> leftrecords = new ArrayList<SAMRecord>();
        Collection<SAMRecord> rightrecords = new ArrayList<SAMRecord>();

        boolean keepgoing = true;
        while (keepgoing) {
            SAMRecord left = nextLeft();
            SAMRecord right = nextRight();
            if (left == null || right == null) {
                break;
            }
            System.err.println("LEFT " + left);
            System.err.println("RIGHT " + right);
            if (left.getReadName().equals(right.getReadName())) {
                leftrecords.add(left);
                rightrecords.add(right);
                System.err.println("MATCH.  Storing.");
            } else {                
                dumpRecords(leftrecords, rightrecords);
                leftrecords.clear();
                rightrecords.clear();                
                leftbuffer.add(left);
                rightbuffer.add(right);
                System.err.println("mismatch.  dumped and cleared\n");
            }
            for (int i = 0; i < leftbuffer.size(); i++) {
                int j = 0;
                while (j < rightbuffer.size() && !leftbuffer.get(i).getReadName().equals(rightbuffer.get(j).getReadName())) {
                    j++;
                }
                if (j == rightbuffer.size()) {
                    continue;
                }
                int k = i;
                int l = j;
                do {
                    leftrecords.add(leftbuffer.get(k++));
                } while (k < leftbuffer.size() && leftbuffer.get(i).getReadName().equals(leftbuffer.get(k).getReadName()));
                do {
                    rightrecords.add(rightbuffer.get(l++));
                } while (l < rightbuffer.size() && rightbuffer.get(i).getReadName().equals(rightbuffer.get(l).getReadName()));
                dumpRecords(leftrecords, rightrecords);
                leftrecords.clear();
                rightrecords.clear();                
                
                while (k-- > 0) {
                    leftbuffer.remove(0);
                }
                while (l-- > 0) {
                    rightbuffer.remove(0);
                }

                i = -1;                
            }
            if (!leftiter.hasNext()) {
                rightbuffer.clear();
            }
            if (!rightiter.hasNext()) {
                leftbuffer.clear();
            }

            System.err.println("li.hn " + leftiter.hasNext() + " lb.size " + leftbuffer.size() + 
                               "ri.hn " + rightiter.hasNext() + " rb.size " + rightbuffer.size());
            keepgoing = (leftiter.hasNext() || leftbuffer.size() > 0) &&
                (rightiter.hasNext() || rightbuffer.size() > 0);
        }
        dumpRecords(leftrecords, rightrecords);
    }
    public static SAMRecord nextLeft() {
        if (leftiter.hasNext()) {
            SAMRecord r = leftiter.next();
            r.setReadName(r.getReadName().replaceAll("/\\d$",""));
            if (r.getReferenceName().equals("*")) {
                return nextLeft();
            } else {
                return r;
            }
        } else if (leftbuffer.size() > 0) {
            return leftbuffer.get(leftbuffer.size() -1);
        } else {
            return null;
        }
    }
    public static SAMRecord nextRight() {
        if (rightiter.hasNext()) {
            SAMRecord r = rightiter.next();
            r.setReadName(r.getReadName().replaceAll("/\\d$",""));
            if (r.getReferenceName().equals("*")) {
                return nextRight();
            } else {
                return r;
            }
        } else if (rightbuffer.size() > 0) {
            return rightbuffer.get(rightbuffer.size() -1);
        } else {
            return null;
        }
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