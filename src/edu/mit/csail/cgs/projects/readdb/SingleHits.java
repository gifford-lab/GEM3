package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;


/** 
 * Represents the list of sorted reads on disk
 */
public class SingleHits extends Hits {
    /**
     * Initializes a Hits object from a file
     */
    public SingleHits (String prefix, int chrom) throws FileNotFoundException, SecurityException, IOException {
        super(chrom,
              getPositionsFname(prefix,chrom),
              getWeightsFname(prefix,chrom), 
              getLaSFname(prefix,chrom));
    }
    public static void writeSingleHits(IntBP positions,
                                       FloatBP weights,
                                       IntBP las,
                                       String prefix,
                                       int chrom) throws IOException {
        String postmp = getPositionsFname(prefix,chrom) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom) + ".tmp";
        RandomAccessFile positionsRAF = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(lastmp,"rw");

        Bits.sendBytes(positions.bb, 0, positions.bb.limit(), positionsRAF.getChannel());
        Bits.sendBytes(weights.bb, 0, weights.bb.limit(), weightsRAF.getChannel());
        Bits.sendBytes(las.bb, 0, las.bb.limit(), lasRAF.getChannel());
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();

        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(getPositionsFname(prefix,chrom));
        (new File(weightstmp)).renameTo(getWeightsFname(prefix,chrom));
        (new File(lastmp)).renameTo(getLaSFname(prefix,chrom));
    }
    public static void writeSingleHits(List<SingleHits> hits,
                                       String prefix, 
                                       int chrom) throws IOException {
        Collections.sort(hits);

        IntBP p = new IntBP(hits.length());
        FloatBP w = new FloatBP(hits.length());
        IntBP l = new IntBP(hits.length());
        for (int i = 0; i < hits.length(); i++) {
            SingleHit h = hits.get(i);
            p.put(i, h.pos);
            w.put(i, h.weight);
            l.put(i, makeLAS(h.length, h.strand));
        }
        writeSingleHits(p,w,l,prefix,chrom);
    }
    private static String getPositionsFname(String prefix, int chrom) {
        prefix + chrom + ".spositions";
    }
    private static String getWeightsFname(String prefix, int chrom) {
        prefix + chrom + ".sweights";
    }
    private static String getLaSFname(String prefix, int chrom) {
        prefix + chrom + ".slas";
    }

}