package edu.mit.csail.cgs.tools.alignments;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.CGSException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.ewok.verbs.io.*;

/**
 * Filters a stream of PSL results to
 *  - require a minimum number of bp of match between query and target
 *  - a minimum percent of matched between query and target
 *  - a maximum number of genomic matches
 *  - only hits that have the maximum number of bp of match
 *  - output only a certain number of matches, selected randomly from the matches meeting the other criteria
 *
 * Usage
 * FilterPSL --species "$SC;SGDv1" --minbpmatch 20 --minpercentbpmatch 80 --maxgenomicmatches 10 --matchestoprint 1 --onlybestscore
 *
 */
public class FilterPSL {
    public static void main(String args[]) throws NotFoundException, IOException, CGSException {
        int minBPMatch = Args.parseInteger(args,"minbpmatch",20);
        double minPercentBPMatch = Args.parseDouble(args,"minpercentbpmatch",0);
        int maxGenomicMatches = Args.parseInteger(args,"maxgenomicmatches",10000);
        int matchesToPrint = Args.parseInteger(args,"matchestoprint",500);
        /* if true, then only accept hits that have the highest score */
        boolean onlyBestScore = Args.parseFlags(args).contains("onlybestscore");
        System.err.println(String.format("Min BP to Match %d, Min %% to match %.3f, Max genomic matches %d, Matches to print %d",
                                         minBPMatch, minPercentBPMatch, maxGenomicMatches, matchesToPrint));

        Pair<Organism,Genome> pair = Args.parseGenome(args);
        Genome genome = pair.cdr();        
        BLATResultFileInputStream stream = new BLATResultFileInputStream(new BufferedReader(new InputStreamReader(System.in)), genome);

        List<PSLHitRegion> buffered = new ArrayList<PSLHitRegion>();
        String lastQID = "";
        
        while (stream.hasNext()) {
            PSLHitRegion r;
            try {
                r = (PSLHitRegion)stream.next();                   
            } catch (Exception e) {
                e.printStackTrace();
                continue;
            }
            if (r == null) {
                continue;
            }


            if (r.getPercentIdentity() < minPercentBPMatch) {
                //                System.err.println("Skipping " + r.getName() + " because " + r.getPercentIdentity() + " < " + minPercentBPMatch);
                continue;
            }
            if (r.getMatches() < minBPMatch) {
                //                System.err.println(String.format("Skipping %s because abs(%d - %d) - %d < %d",
                                                 //                                                 r.getName(), r.getQueryEnd(), r.getQueryStart(), r.getNumMismatches(), minBPMatch));
                continue;
            }
            if (!r.getName().equals(lastQID)) {
                //                System.err.println("Considering " + lastQID + " at " + buffered.size());
                if (buffered.size() <= maxGenomicMatches) {
                    printHits(buffered,
                              matchesToPrint,
                              onlyBestScore);
                }

                buffered.clear();
            }
            buffered.add(r);
            lastQID = r.getName();
        }        
        if (buffered.size() <= maxGenomicMatches) {
            printHits(buffered,
                      matchesToPrint,
                      onlyBestScore);
        }
    }

    public static void printHits (List<PSLHitRegion> hits,
                                  int tokeep,
                                  boolean onlybestscore) {
        int score = 0;
        if (onlybestscore) {
            for (PSLHitRegion r : hits) {
                if (r.getMatches() > score) {
                    score = r.getMatches();
                }
            }
            int i = 0;
            while (i < hits.size()) {
                if (hits.get(i).getMatches() < score) {
                    hits.remove(i);
                } else {
                    i++;
                }
            }
        }

        while (hits.size() > tokeep) {
            int toRemove = (int)(Math.floor(Math.random() * hits.size()));
            hits.remove(toRemove);
        }
        for (PSLHitRegion h : hits) {
            System.out.println(h.toString());
        }
    }
}