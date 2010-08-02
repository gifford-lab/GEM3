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
 * java RegionListOverlap --species "$MM;mm9" --one fileone.txt --two filetwo.txt --colone 0 --coltwo 3 --window 50
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

        readFile(genome, fone, colone, one);
        readFile(genome, ftwo, coltwo, two);

        for (String chrom : one.keySet()) {
            if (!two.containsKey(chrom)) { continue;}
            List<Region> lone = one.get(chrom);
            List<Region> ltwo = two.get(chrom);
            Collections.sort(lone);
            Collections.sort(ltwo);
            for (Region orig : lone) {
                Region r = orig.expand(window,window);
                for (Region o : ltwo) {
                    if (r.overlaps(o)) {
                        System.out.println(orig.toString());
                        break;
                    }
                }
            }
        }
    }
    public static void readFile(Genome genome, String fname, int column, Map<String,List<Region>> map) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(fname));
        String line = null;
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("\\t");
            Region r = Region.fromString(genome, pieces[column]);
            if (!map.containsKey(r.getChrom())) {
                map.put(r.getChrom(), new ArrayList<Region>());
            }
            map.get(r.getChrom()).add(r);
        }
        reader.close();
    }

}