package edu.mit.csail.cgs.tools.regions;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Reads regions from two tab delimited files (you specify 
 * the columns, default is first column) and determines
 * how many from the first file are present in the second
 * given some mismatch window.
 *
 * java RegionListOverlap --species "$MM;mm9" --one fileone.txt --two filetwo.txt --colone 0 --coltwo 3 --window 50 [--stats]
 *
 * --stats says to use --one as a test set and --two as a gold standard set to report the TP and FP rates for one
 */

public class RegionListOverlap {

    public static void main(String args[]) throws Exception {
        Map<String,List<Region>> one = new HashMap<String,List<Region>>();
        Map<String,List<Region>> two = new HashMap<String,List<Region>>();
        String fone = Args.parseString(args,"one",null);
        String ftwo = Args.parseString(args,"two",null);
        int colone = Args.parseInteger(args,"colone",0);
        int coltwo = Args.parseInteger(args,"coltwo",0);
        int window = Args.parseInteger(args,"window",0);
        Genome genome = Args.parseGenome(args).cdr();
        boolean stats = Args.parseFlags(args).contains("stats");

        readFile(genome, fone, colone, one);
        readFile(genome, ftwo, coltwo, two);

        int overlap = 0;
        int onecount = count(one);
        int twocount = count(two);

        for (String chrom : one.keySet()) {
            if (!two.containsKey(chrom)) { continue;}
            List<Region> lone = one.get(chrom);
            List<Region> ltwo = two.get(chrom);
            Collections.sort(lone);
            Collections.sort(ltwo);
            for (Region orig : lone) {
                Region r = orig.expand(window,window);
                boolean found = false;
                for (Region o : ltwo) {
                    if (r.overlaps(o)) {
                        if (stats) {
                            System.out.println("TP\t" + orig);
                        } else {
                            System.out.println(orig.toString());
                        }
                        overlap++;
                        found = true;
                        break;
                    }
                }
                if (stats && !found) {
                    System.out.println("FP\t" + orig);
                }

            }
        }
        if (stats) {
            int tp = overlap;
            int fp = onecount - overlap;
            int fn = twocount - overlap;
            double tprate = ((double)tp) / onecount;
            double fnrate = ((double)overlap) / twocount;
            System.out.println(String.format("n=%d  tp=%d  fp=%d  fn=%d   tpr=%.2f   fnr=%.2f",
                                             onecount, tp,fp,fn,tprate,fnrate));
        }
    }
    public static int count(Map<String,List<Region>> map) {
        int count = 0;
        for (String chrom : map.keySet()) {
            count += map.get(chrom).size();
        }
        return count;
    }
    public static void readFile(Genome genome, String fname, int column, Map<String,List<Region>> map) throws IOException {
        if (map == null) {
            throw new NullPointerException("Can't give me a null map");
        }
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = null;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");
            Region r = Region.fromString(genome, pieces[column]);
            if (r == null) {
                System.err.println("Couldn't parse " + pieces[column] + " in " + fname);
                continue;
            }
            if (!map.containsKey(r.getChrom())) {
                map.put(r.getChrom(), new ArrayList<Region>());
            }
            map.get(r.getChrom()).add(r);
        }
        reader.close();
    }

}