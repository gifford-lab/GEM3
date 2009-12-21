package edu.mit.csail.cgs.datasets.fragdist;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import cern.colt.matrix.*;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.regression.Regression;

/** 
 * TimeToLen processes Agilent Bioanalyzer output.  It uses
 * the ladder file (which maps time to a fragment size in bp)
 * to map files of observations from the time domain to the
 * length domain.
 *
 * Usage:
 *  java edu.mit.csail.cgs.datasets.fragdist --ladder ladder.time --minlength 50 --maxlength 1400 --sample lane1.csv --sample lane2.csv
 * will produce lane1.raw and lane2.raw.  The ladder.time file comes from PeakID.  The resulting lane.raw files
 * can be processed with ConToInt.
 *
 */
public class TimeToLen {    

    public static void main(String args[]) throws Exception {
        String ladderfname = null;
        int minlength = 0;
        int maxlength = 100000;
        ArrayList<String> samplefnames = new ArrayList<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--ladder")) {
                ladderfname = args[++i];
            }
            if (args[i].equals("--sample")) {
                samplefnames.add(args[++i]);
            }
            if (args[i].equals("--minlength")) {
                minlength = (int)Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--maxlength")) {
                maxlength = (int) Integer.parseInt(args[++i]);
            }


        }
        TimeToLen object = new TimeToLen(ladderfname,
                                         samplefnames);
        object.runAndPrint(minlength,maxlength);
    }

    private int ladderMaxSize, ladderMinSize;
    private double ladderSlope, ladderIntercept, ladderMinTime, ladderMaxTime;
    private List<Double> ladderTimes;
    private List<Integer> ladderLengths;
    private List<String> samplefnames;
       
    public TimeToLen(String ladderfname,
                     List<String> samplefnames) throws IOException {
        readLadder(ladderfname);
        this.samplefnames = samplefnames;
    }
    public void readLadder(String ladderfname) throws IOException, FileNotFoundException {
        ArrayList<Double> times = new ArrayList<Double>();
        ArrayList<Integer> lengths = new ArrayList<Integer>();
        BufferedReader reader = new BufferedReader(new FileReader(ladderfname));
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.matches("\\d+\\t[\\d\\.]+")) {
                String pieces[] = line.split("\\t");                
                lengths.add(Integer.parseInt(pieces[0]));
                times.add(Double.parseDouble(pieces[1]));
            }
        }
        reader.close();
        ladderMinSize = lengths.get(0) - 1;   // don't want to try to take log(0) later on
        ladderMaxSize = lengths.get(lengths.size() - 1);
        ladderMinTime = times.get(0);
        ladderMaxTime = times.get(times.size() - 1);

        DoubleMatrix2D dep = DoubleFactory2D.dense.make(lengths.size(),1);
        DoubleMatrix2D indep = DoubleFactory2D.dense.make(lengths.size(),2);
        for (int i = 0; i < lengths.size(); i++) {
            dep.set(i,1,Math.log(lengths.get(i) - ladderMinSize));
            indep.set(i,0,1);
            indep.set(i,1,times.get(i) - ladderMinTime);
        }
        DoubleMatrix2D coeffs = Regression.linear(indep,dep);
        ladderIntercept = coeffs.get(0,0);
        ladderSlope = coeffs.get(1,0);

        DoubleMatrix2D pred = Regression.predict(indep,coeffs);
        Pair<Double,Double> score = Regression.score(pred,dep);
        System.err.println("LadderIntercept= " + ladderIntercept + "   ladderSlope=" + ladderSlope);
        System.err.println("Ladder r2=" + score.getLast());
        for (int i = 0; i < lengths.size(); i++) {
            System.err.println("  time=" + times.get(i) + "  length=" + Math.exp(dep.get(i,0)) + "  predlen=" + Math.exp(pred.get(i,0)));
        }


        ladderLengths = lengths;
        ladderTimes = times;
    }        

//     public Map<Double,Double> timeToLen(Map<Double,Double> times) {
//         HashMap<Double,Double> out = new HashMap<Double,Double>();
//         int ladderindex = 0;
//         List<Double> keys = new ArrayList<Double>();
//         keys.addAll(times.keySet());
//         Collections.sort(keys);
//         double k = 0;
//         k = Math.log(ladderLengths.get(1) - ladderLengths.get(0)) / (ladderTimes.get(1) - ladderTimes.get(0));
//         for (Double d : keys) {
//             if (d > ladderTimes.get(ladderindex+1)) {
//                 if (ladderindex + 2 < ladderTimes.size()) {
//                     ladderindex++;
//                     k = Math.log(ladderLengths.get(ladderindex+1) - ladderLengths.get(ladderindex)) / 
//                         (ladderTimes.get(ladderindex+1) - ladderTimes.get(ladderindex));                    
//                 } else {
//                     break;
//                 }
//             }
//             double length = ladderLengths.get(ladderindex) + Math.exp(k * (d - ladderTimes.get(ladderindex)));
//             out.put(length,times.get(d));
//         }
//         return out;        
//     }

    public Map<Double,Double> timeToLen(Map<Double,Double> times) {
        HashMap<Double,Double> out = new HashMap<Double,Double>();
        int ladderindex = 0;
        List<Double> keys = new ArrayList<Double>();
        keys.addAll(times.keySet());
        Collections.sort(keys);
        for (Double d : keys) {
            if (d > ladderMaxTime) {
                break;
            }
            out.put(ladderMinSize + Math.exp(ladderIntercept + ladderSlope * (d - ladderMinTime)),times.get(d));
        }
        return out;        
    }
    

    public void runAndPrint(int minlength, int maxlength) throws IOException, FileNotFoundException{
        for (String fname : samplefnames) {
            String outfname = fname.replaceAll("\\.\\w{2,4}$",".raw");
            if (outfname.equals(fname)) {
                throw new RuntimeException("NAMES ARE SAME " + fname);
            }
            System.err.println("Writing to " + outfname);
            PrintWriter writer = new PrintWriter(outfname);
            Map<Double,Double> times = CSV.readFile(fname);
            Map<Double,Double> lengths = timeToLen(times);
            List<Double> keys = new ArrayList<Double>();
            keys.addAll(lengths.keySet());
            Collections.sort(keys);
            for (Double d : keys) {
                if (d < minlength || d > maxlength) {
                    continue;
                }
                writer.println(d + "\t" + lengths.get(d));
            }
            writer.close();
            
        }
    }

    

}
