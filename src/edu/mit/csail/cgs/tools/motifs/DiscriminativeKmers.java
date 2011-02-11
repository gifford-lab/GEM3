package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.*;


public class DiscriminativeKmers {

    private Genome genome;
    private double minfoldchange;
    private int k, mask;
    private Binomial binomial;
    private SequenceGenerator seqgen;
    private int randombgcount = 1000,  // number of random background regions to pic
        randombgsize = 100, // size of random background regions
        parsedregionexpand; // expand input regions by this much on either side
    private Map<String,char[]> foreground, background;  // foreground and background sequences

    public static long addChar(long existing, long mask, char newchar) {
        existing <<= 2;
        switch (newchar) {
        case 'A':
        case 'a':
            existing += 0;
            break;
        case 'C':
        case 'c':
            existing += 1;
            break;
        case 'G':
        case 'g':
            existing += 2;
            break;
        case 'T':
        case 't':
            existing += 3;
            break;

        default:
            break;
        }
        return existing & mask;
    }
    public static String longToString(long l, int len) {
        StringBuilder builder = new StringBuilder();
        while (len-- > 0) {
            int byte = l & 3;
            l >>= 2;
            switch (byte) {
            case 0:
                builder.append("A");
                break;
            case 1:
                builder.append("C");
                break;
            case 2:
                builder.append("G");
                break;
            case 3:
                builder.append("T");
                break;
            default:
                break;
            }
        }
        builder.reverse();
        return builder.toString();
    }
    public static Map<Long,Integer> count(char[] chars, 
                                          int k,
                                          long mask,
                                          Map<Long,Integer> map) {
        if (map == null) {
            map = new HashMap<Long,Integer>();
        }
        long l = 0;
        for (int i = 0; i < k-1; i++) {
            l = addChar(l,chars[i]);
        }
        while (int i = k; i < chars.length; i++) {
            l = addChar(l,mask,chars[i]);
            if (!map.containsKey(l)) {
                map.put(l,1);
            } else {
                map.put(l,map.get(l)+1);
            }            
        }
        return map;
    }
    public DiscriminativeKmers () {
        seqgen = new SequenceGenerator();
    }
    public setK(int k) {
        this.k = k;
        mask = 0;
        for (int i = 0; i < k; i++) {
            mask <<= 1;
            mask = mask | 1;
        }
    }
    public void parseArgs(String args[]) {
        int k = Args.parseInteger(args,"k",10);
        setK(k);
        minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        parseregionexpand = Args.parseInteger(args,"expand",30);
        randombgcount = Args.parseInteger(args,"randombgcount",1000);
        randombgsize = Args.parseInteger(args,"randombgsize",100);
        genome = Args.parseGenome(args).cdr();
        String firstfname = null, secondfname = null;
        firstfname = Args.parseString(args,"first",null);
        secondfname = Args.parseString(args,"second",null);
        if (firstfname == null) {
            System.err.println("No --first specified.  Reading from stdin");
            foreground = CompareEnrichment.readRegions(genome,
                                                       new BufferedReader(new InputStreamReader(System.in)), 
                                                       parsedregionexpand, 
                                                       matchedbg ? background : null);
        } else {
            if (firstfname.matches(".*\\.fasta") ||
                firstfname.matches(".*\\.fa")) {
                foreground = readFasta(new BufferedReader(new FileReader(firstfname)));
            } else {
                foreground = CompareEnrichment.readRegions(genome,
                                                           new BufferedReader(new FileReader(firstfname)), 
                                                           parsedregionexpand,
                                                           matchedbg ? background : null);
            }
        }
        if (secondfname == null) {
            System.err.println("No background file given.  Generating " + randombgcount + " regions of size " + randombgsize);
            background = CompareEnrichment.randomRegions(genome,
                                                         randombgcount,
                                                         randombgsize);
        } else {
            if (secondfname.matches(".*\\.fasta") ||
                secondfname.matches(".*\\.fa")){
                background = CompareEnrichment.readFasta(new BufferedReader(new FileReader(secondfname)));
            } else {
                background = CompareEnrichment.readRegions(genome, new BufferedReader(new FileReader(secondfname)), parsedregionexpand,null);
            }
        }
    }
    public void run() {
        Map<Long,Integer> fgcounts = new HashMap<Long,Integer>();
        Map<Long,Integer> bgcounts = new HashMap<Long,Integer>();
        for (char[] chars : foreground.values()) {
            count(chars, k, mask, fgcounts);
        }
        for (char[] chars : background.values()) {
            count(chars, k, mask, bgcounts);
        }
        int fgsum = 0;
        int bgsum = 0;
        for (long l : fgcounts.keySet()) {
            fgsum += fgcounts.get(l);
        }
        for (long l : bgcounts.keySet()) {
            bgsum += bgcounts.get(l);
        }

        for (long l : fgcounts.keySet()) {
            int bgcount = (bgcounts.contains(l) ? bgcounts.get(l) : 0) + 1;
            double bgprob = ((double)bgcount) / ((double)bgsum);
            int fgcount = fgcounts.get(l);
            double fgprob = ((double)fgcount) / ((double)fgsum);

            if (fgprob > bgprob * minfoldchange) {
                
            }

        }



    }


}