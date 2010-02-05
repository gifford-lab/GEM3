package edu.mit.csail.cgs.tools.alignments;

import java.io.*;
import edu.mit.csail.cgs.utils.parsing.alignment.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Reads SAM formatted alignment output on STDIN and produces
 * a list of the SNP, insertion, and deletion rate at
 * each position in the read
 *
 * Usage:
 * SAM2ErrorRate --n 36 < alignment.sam
 */

public class SAM2ErrorRate extends SAMHandler {

    private int n, skipped, noalign;
    // the counts of what we saw at each position
    private int[] correct, snp, insertion, deletion;
    // represents current read
    private boolean unchanged[], skipThisRead;
    private int stage;

    public static void main(String[] args) throws IOException {
        SAM2ErrorRate rate = new SAM2ErrorRate(new BufferedReader(new InputStreamReader(System.in)));
        rate.parseArgs(args);
        rate.parse();
        rate.report();
    }

    public SAM2ErrorRate(BufferedReader r) {
        super(r);
    }

    public void parseArgs(String[] args) {
        n = Args.parseInteger(args,"n",-1);
        if (n == -1) {
            throw new IllegalArgumentException("Must supply --n as the read length");
        }
        correct = new int[n];
        snp = new int[n];
        insertion = new int[n];
        deletion = new int[n];
        unchanged = new boolean[n];
        skipped = 0;
        noalign = 0;
        for (int i = 0; i < n; i++) {
            correct[i] = 0;
            snp[i] = 0;
            insertion[i] = 0;
            deletion[i] = 0;
        }
    }

    public void handleRecord(SAMRecord record) {
        if ((record.flags & SAMParser.UNMAPPEDQUERY) != 0) {
            noalign++;
            return;
        }

        for (int i = 0; i < n; i++) {
            unchanged[i] = true;
        }
        skipThisRead = false;
        stage = 0;
        for (int i = 0; i < record.nfields; i++) {
            if (record.tags[i].equals("MD")) {                    
                parseMDTag(record.values[i]);
            }
        }
        parseCIGAR(record.cigar);       
        stage = 1;
        if (skipThisRead) {
            skipped++;
            return;
        }
        parseCIGAR(record.cigar);              
    }
    public void handleMDLetters(int readpos, String letters) {
        for (int i = 0; i < letters.length(); i++) {
            unchanged[i+readpos] = false;
        }
    }
    public void handleCIGARPart(char type, int count, int readpos) {
        if (stage == 0) {
            if (type == 'I') {
                for (int i = 0; i < count; i++) {
                    for (int j = unchanged.length-2; j >= readpos; j--) {
                        unchanged[j+1] = unchanged[j];
                    }
                    unchanged[readpos] = false;
                }        
            }    
            if (type == 'S' || type == 'N' || type == 'H' || type == 'P') {
                skipThisRead = true;
            }
            return;
        } else {
            if (type == 'M') {
                for (int j = 0; j < count; j++) {
                    if (unchanged[readpos+j]) {
                        correct[readpos+j]++;
                    } else {
                        snp[readpos+j]++;
                    }
                }
            }
            if (type == 'D') {
                deletion[readpos]++;
            }
            if (type == 'I') {
                insertion[readpos]++;
            }
        }
    }
    
    public void report() {
        System.err.println("No Alignment for " + noalign + " and skipped " + skipped);
       for (int i = 0; i < n; i++) {
            int total = correct[i] + snp[i] + insertion[i] + deletion[i];
            System.out.println(String.format("%d : %d total, %d correct, %d snp, %d ins, %d del",
                                             i,total,correct[i], snp[i], insertion[i], deletion[i]));
        }

    }
    

}