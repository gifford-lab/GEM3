package edu.mit.csail.cgs.datasets.fragdist;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/** Converts the ladder.csv file to a ladder.time file by
 * identifying a set of known peak locations.
 *
 * Usage:
 * java edu.mit.csail.psrg.Datasets.fragdist.PeakID < Ladder.csv > ladder.time
 *
 * The default is to assume the 1500bp bioanalyzer chip.  --long means
 * to use the marker locations for the 15kb chip.
 */
public class PeakID {

    public static final int[] ladderPeaksShort = {15,25,50,100,150,200,300,400,500,700,850,1000,1500};
    public static final int[] ladderPeaksLong = {50,100,300,500,700,1000,1500,2000,3000,5000,7000,10380};
    
    public static void main(String args[]) throws Exception {
        boolean isShort = true, resultsfile = false;
        int[] levels;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--long")) {
                isShort = false;
            }
            if (args[i].equals("--results")) {
                resultsfile = true;
            }

        }
        if (isShort) {
            levels = ladderPeaksShort;
        } else {
            levels = ladderPeaksLong;
        }

        List<Double> peaks = new ArrayList<Double>(), peakvals = new ArrayList<Double>();
        Map<Double,Double> times;
        if (resultsfile) {
            readResultsForLadder(new BufferedReader(new InputStreamReader(System.in)),
                                 peakvals,
                                 peaks);
        } else {
            times = CSV.readFile(new BufferedReader(new InputStreamReader(System.in)));
            ArrayList<Double> keys = new ArrayList<Double>();
            keys.addAll(times.keySet());
            Collections.sort(keys);
            double prev = 0, peak = 0, peakval = 0; int descend = 0;
            boolean ladderpeak = false;        
            for (double i : keys) {
                if (times.get(i) > 0) {
                    if (times.get(i) > prev) {
                        descend = 0;
                        peak = i;
                        peakval = times.get(i);
                        if (peakval > 10) {
                            ladderpeak = true;
                        }                    
                    } else {
                        descend++;
                        if (descend > 3 && ladderpeak) {
                            peaks.add(peak);
                            peakvals.add(peakval);
                            ladderpeak = false;
                        }
                        
                    }
                    
                }
                prev = times.get(i);
            }
        }
        if (peaks.size() == levels.length) {
            for (int l = 0; l < peaks.size(); l++) {
                System.out.println(levels[l] + "\t" + peaks.get(l));
            }
        } else {
            System.err.println("WARNING: can't find correct number of peaks");
            System.err.println("Expected " + levels.length + " but found " + peaks.size());
            for (int l = 0; l < Math.min(peaks.size(), levels.length); l++) {
                System.out.println(levels[l] + "\t" + peaks.get(l));
            }

        }
    }

    public static void readResultsForLadder(BufferedReader reader, List<Double> distances, List<Double> times) throws IOException {
        String line;
        boolean inLadder = false, inTable = false;
        while ((line = reader.readLine()) != null) {
            if (line.matches("Sample Name,Ladder.*")) {
                inLadder = true;
            } else if (inLadder && line.matches("Peak Table.*")) {
                inTable = true;
            } else if (inTable) {
                if (line.matches("Size.*")) {
                    continue;
                }
                String pieces[] = line.split(",");
                if (pieces.length <= 1) {
                    inTable = false;
                    break;
                }

                pieces[0] =pieces[0].replaceAll("\\D","");
                Double length = Double.parseDouble(pieces[0]);
                Double time = Double.parseDouble(pieces[5]);
                distances.add(length);
                times.add(time);
            }
        }
    }
}
