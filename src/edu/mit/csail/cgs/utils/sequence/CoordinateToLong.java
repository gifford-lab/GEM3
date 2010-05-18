package edu.mit.csail.cgs.utils.sequence;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Point;
import java.util.*;

/**
 * Maps coordinates in a genome to ints and back again.
 */

public class CoordinateToLong {

    private Genome genome;
    private List<String> chroms;
    private Map<String,Long> offsets;
    private long[] sorted;

    public CoordinateToLong(Genome g) {
        genome = g;
        chroms = new ArrayList<String>();
        offsets = new HashMap<String,Long>();
        Map<String,Integer> lengths = genome.getChromLengthMap();

        chroms.addAll(lengths.keySet());
        Collections.sort(chroms);
        long offset = 0;
        for (String s : chroms) {
            offsets.put(s,offset);
            offset += lengths.get(s);
        }
        sorted = new long[chroms.size()];
        for (int i = 0; i < sorted.length; i++) {
            sorted[i] = offsets.get(chroms.get(i));
        }
    }
    public long convert(String c, int pos) {
        return offsets.get(c) + pos;
    }
    public long convert(Point p) {
        return convert(p.getChrom(), p.getLocation());
    }
    public Point convert(long l) {
        int i = Arrays.binarySearch(sorted, l);
        if (i < 0) {
            i = (i+1)*-1 - 1;
        }
        String s = chroms.get(i);
        return new Point(genome, s, (int)(l - sorted[i]));
    }
    public long getOffset(long l) {
        int i = Arrays.binarySearch(sorted, l);
        if (i < 0) {
            i = (i+1)*-1 - 1;
        }
        return l - sorted[i];
    }
    public String getChrom(long l) {
        int i = Arrays.binarySearch(sorted, l);
        if (i < 0) {
            i = (i+1)*-1 - 1;
        }
        String s = chroms.get(i);
        return s;
    }


}