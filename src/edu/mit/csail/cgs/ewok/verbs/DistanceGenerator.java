package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;

import java.util.*;

/**
 * DistanceGenerator takes two Expander<Region,Region>, a source and a
 * sink.  When execute() is presented with a Region, DistanceGenerator
 * applies the two expanders to that region to generate a set of
 * sources and a set of sinks.  It then tries to map each source to
 * the nearest sink without reusing any sink.  There is a maximum
 * distance (distanceLimit) for the matching.
 *
 * execute() doesn't allow you to determine the distance for any given
 * source-sink; instead, it returns an iterator over all the distances
 * with negative values indicating a source that couldn't be mapped
 * (no sinks left) or was too far.
 *
 * This class is useful, for example, for mapping binding events to
 * motifs.
 */
public class DistanceGenerator implements Expander<Region,Integer> {

    private Expander<Region,Region> sourcemapper, sinkmapper;
    /* how far is "toofar" */
    private int distanceLimit = 1000;
    /* how much should we expand the input Region
       when presenting it to the sinkmapper */
    private int expandInputRegion;

    public DistanceGenerator(Expander<Region,Region> sourcemapper,
                             Expander<Region,Region> sinkmapper) {
        this.sourcemapper = sourcemapper;
        this.sinkmapper = sinkmapper;
        expandInputRegion = 1000;
        distanceLimit = 1000;
    }
    public DistanceGenerator(Expander<Region,Region> sourcemapper,
                             Expander<Region,Region> sinkmapper,
                             int expandInputRegion,
                             int distanceLimit) {
        this.sourcemapper = sourcemapper;
        this.sinkmapper = sinkmapper;
        this.expandInputRegion = expandInputRegion;
        this.distanceLimit = distanceLimit;
    }

    

    /** 
     * returns a distance of -1 for any source that can't be
     * matched to a sink and -2 for anything that's too far 
     */
    public Iterator<Integer> execute(Region region) {
        ArrayList<Integer> sourcelist, sinklist, distances;
        Integer[] sources, sinks;
        sourcelist = new ArrayList<Integer>();
        sinklist = new ArrayList<Integer>();
        Iterator<Region> sourceIter = sourcemapper.execute(region);
        Region sinkRegion = region.expand(-1*expandInputRegion,expandInputRegion);
        Iterator<Region> sinkIter = sinkmapper.execute(sinkRegion);
        while (sourceIter.hasNext()) {
            Region r = sourceIter.next();
            sourcelist.add((r.getStart() + r.getEnd())/2);
        }
        while (sinkIter.hasNext()) {
            Region r = sinkIter.next();
            sinklist.add((r.getStart() + r.getEnd())/2);
        }
//         System.err.println("Regin is " + region);
//         System.err.println("Sources are " + sourcelist);
//         System.err.println("Sinks are " + sinklist);
        sources = sourcelist.toArray(new Integer[0]);
        sinks = sinklist.toArray(new Integer[0]);
        Arrays.sort(sources);
        Arrays.sort(sinks);
        distances = new ArrayList<Integer>();
        int todo = sources.length;
        int snkleft = sinks.length;
        for (int i = 0; i < sources.length - 1; i++) {
            if (sources[i] == sources[i+1]) {
                sources[i] = -1;
                todo--;
            }
        }
        int before = 0, after = 0;
        while ((todo > 0) && (snkleft > 0)) {
            boolean anything = false;
            for (int i = 0; i < sources.length; i++) {
                if (sources[i] == -1) {continue;}
                while ((after <= sinks.length - 1) &&
                       ((sinks[after] < sources[i]) ||
                        (sinks[after] == -1)))  {after++;}
                before = after - 1;
                while ((before > 0) &&
                       (sinks[before] == -1)) {before--;}                                
                int d1 = Integer.MAX_VALUE, d2 = Integer.MAX_VALUE, d3 = Integer.MAX_VALUE, d4 = Integer.MAX_VALUE;
                int prevsrc = i - 1, nextsrc = i + 1;
                while ((prevsrc > 0) &&
                       (sources[prevsrc] == -1)) {prevsrc--;}
                while ((nextsrc < sources.length - 1) &&
                       (sources[nextsrc] == -1)) {nextsrc++;}
                if ((before >= 0) &&
                    (prevsrc >= 0) &&
                    (prevsrc != i) &&
                    (sources[prevsrc] != -1) &&
                    (sinks[before] != -1)) {
                    d1 = Math.abs(sources[prevsrc] - sinks[before]);
                }
                if ((before >= 0) &&
                    (sinks[before] != -1)) {
                    d2 = Math.abs(sources[i] - sinks[before]);
                }
                if ((after < sinks.length) &&
                    (sinks[after] != -1)) {
                    d3 = Math.abs(sources[i] - sinks[after]);
                }
                if ((after < sinks.length) &&
                    (sinks[after] != -1) &&
                    (nextsrc < sources.length) &&
                    (nextsrc != i) &&
                    (sources[nextsrc] != -1)) {
                    d4 = Math.abs(sources[nextsrc] - sinks[after]);
                }
                //                System.err.println("sources[i] = " + sources[i] + "  sinks are " + Arrays.toString(sinks));
                //                System.err.println("d1=" + d1 + "  d2=" + d2 + "  d3="+ d3 + "  d4=" + d4);
                if (d2 < d3) {
                    if (d2 < d1) {
                        if (d2 < distanceLimit) {
                            //                            System.out.println(d2 + "\t" + region.getChrom() + "\t" + sources[i] + "\t" + sinks[before]);
                            distances.add(d2);
                        } else {
                            distances.add(-2);
                        }
                        sources[i] = -1;
                        sinks[before] = -1;
                        todo--;
                        snkleft--;                            
                        anything = true;
                    }
                } else {
                    if (d3 < d4) {
                        if (d3 < distanceLimit) {
                            distances.add(d3);
                            //                            System.out.println(d3 + "\t" + region.getChrom() + "\t" + sources[i] + "\t" + sinks[after]);
                        } else {
                            distances.add(-2);
                        } 
                        sources[i] = -1;
                        sinks[after] = -1;
                        todo--;
                        snkleft--;
                        anything = true;
                    }
                }
            }
            while (todo > 0) {
                distances.add(-1);
                todo--;
            }
            if (!anything) {
                System.err.println("Didn't do anything");
                for (int i = 0; i < sources.length; i++) {
                    System.err.print(sources[i] + " ");
                }
                System.err.println("");
                for (int i = 0; i < sinks.length; i++) {
                    System.err.print(sinks[i] + " ");
                }
                System.err.println("");
                throw new RuntimeException("Stuck!");
            }
        }
        //        System.err.println("DG is returning " + distances.size());
        return distances.iterator();
    }


}
