package edu.mit.csail.cgs.utils.io.parsing.alignment;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * Superclass/driver for event-based parsing of SAM alignment output
 */
public class SAMHandler {

    public static final String VERSION = "$Id$";

    protected BufferedReader reader;
    protected SAMParser parser;
    protected SAMRecord record;

    /* the current chromosome name (or target sequence name) */
    protected String currentChrom;
    /* the current position */
    protected int currentPos;

    /* parse optional fields from the input? */
    protected boolean wantOptionalFields;

    private Pattern cigarpattern, initialpattern, mdpattern;

    public SAMHandler(BufferedReader r) {
        parser = new SAMParser();
        reader = r;
        record = new SAMRecord();
        wantOptionalFields = true;
        cigarpattern = Pattern.compile("(\\d+[MIDNSHP])");
        initialpattern = Pattern.compile("^(\\d+).*");
        mdpattern = Pattern.compile("(\\^?[ACGTN]+\\d+)");

    }

    /* ------------------------------------------------------------
       These functions are for the main parsing loop
    */
    /**
     * handle the transition from parsing results on one chromosome
     * to those on a new chromosome.
     */
    public void newChrom(String chromName) {}
    /**
     * handle the transition from parsing results starting at one
     * genome coordinate to a new coordinate.
     */
    public void newPos(int pos) {}
    /** 
     * called after reading a record from the input but before doing
     * anything else.  Return false if you want to ignore the record
     */
    public boolean shouldHandleRecord(SAMRecord record) {return true;}
    
    public void parse() throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            parseLine(line);
        }
    }
    public void parseLine(String line) {
        try {
            parser.parseFromLine(line, record, true);
        } catch (Exception e) {
            System.err.println(e.toString() + " from " + line);
            e.printStackTrace();
            return;
        }
        if (!shouldHandleRecord(record)) { return;}
        if (!record.tname.equals(currentChrom)) {
            newChrom(record.tname);
            currentChrom = record.tname;
            currentPos = -100;       
        }
        if (record.pos != currentPos) {
            newPos(record.pos);
            currentPos = record.pos;
        }
        handleRecord(record);
    }
    /**
     * Do whatever we want to do for the current record.  Might call parseCIGAR or parseMD
     */
    public void handleRecord(SAMRecord record) {}

    /* -----------------------------------------------------------
       These functions are for parsing the parts of the record
    */

    public void parseCIGAR(String cigar) {
        Matcher m;        
        m = cigarpattern.matcher(cigar);
        int readpos = 0;
        while (m.find()) {
            String group = m.group();
            if (group == null) { continue;}
            char type = group.charAt(group.length() - 1);
            int count = Integer.parseInt(group.substring(0,group.length() - 1));
            handleCIGARPart(type, count, readpos);
            if (type == 'M' || type == 'I') {
                readpos += count;
            }        
        }
    }
    /**
     * handle part of the CIGAR string.  type (M,I,D) and count are from the string.
     * readpos is the zero-based position at which the letters described by
     * type and count start in the read.
     */
    public void handleCIGARPart(char type, 
                                int count,
                                int readpos) {}


    public void parseMDTag (String tag) {
        int readpos = 0;
        Matcher m = initialpattern.matcher(tag);
        if (!m.find()) {
            System.err.println("Couldn't match initial pattern against " + tag);
            return;
        }
        String initialgroup = m.group(1);
        readpos = Integer.parseInt(initialgroup);
        handleMDSame(0,readpos);
        m = mdpattern.matcher(tag.substring(initialgroup.length()));
        while (m.find()) {
            String group = m.group();
            if (group == null) { continue;}
            if (group.charAt(0) == '^') {
                int k;
                for (k = 1; k < group.length(); k++) {
                    char c = group.charAt(k);
                    if (!Character.isLetter(c)) {break;}
                }         
                handleMDDeletion(readpos, group.substring(1, k));                       
//                 System.err.println(String.format("Saw deletion of %s at %d in %s, %s, %s\n",
//                                                  group.substring(1,k), readpos, record.seq, record.cigar, tag));
                int same = Integer.parseInt(group.substring(k));
                handleMDSame(readpos, same);
                readpos += same;
            } else {
                int k;
                for (k = 0; k < group.length(); k++) {
                    char c = group.charAt(k);
                    if (!Character.isLetter(c)) {break;}
                }
                handleMDLetters(readpos, group.substring(0,k));
                readpos += k;
                int same = Integer.parseInt(group.substring(k));
                handleMDSame(readpos, same);
                readpos += same;
            }
        }
    }
    public void handleMDDeletion(int readpos, String deletedBases) {}
    public void handleMDSame(int readpos, int count) {}
    public void handleMDLetters(int readpos, String letters) {}

}


 