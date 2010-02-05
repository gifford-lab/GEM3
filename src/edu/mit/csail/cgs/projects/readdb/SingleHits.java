package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.List;
import java.util.Collections;

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

        if (positions.bb.order() != ByteOrder.nativeOrder()) {
            Bits.flipByteOrder(positions.ib);            
            positions.bb.order(ByteOrder.nativeOrder());
        }
        if (weights.bb.order() != ByteOrder.nativeOrder()) {
            Bits.flipByteOrder(weights.fb);
            weights.bb.order(ByteOrder.nativeOrder());
        }
        if (las.bb.order() != ByteOrder.nativeOrder()) {
            Bits.flipByteOrder(las.ib);
            las.bb.order(ByteOrder.nativeOrder());
        }

        Bits.sendBytes(positions.bb, 0, positions.bb.limit(), positionsRAF.getChannel());
        Bits.sendBytes(weights.bb, 0, weights.bb.limit(), weightsRAF.getChannel());
        Bits.sendBytes(las.bb, 0, las.bb.limit(), lasRAF.getChannel());
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();

        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom)));
    }
    public static void writeSingleHits(List<SingleHit> hits,
                                       String prefix, 
                                       int chrom) throws IOException {
        Collections.sort(hits);

        IntBP p = new IntBP(hits.size());
        FloatBP w = new FloatBP(hits.size());
        IntBP l = new IntBP(hits.size());
        for (int i = 0; i < hits.size(); i++) {
            SingleHit h = hits.get(i);
            p.put(i, h.pos);
            w.put(i, h.weight);
            l.put(i, makeLAS(h.length, h.strand));
        }
        writeSingleHits(p,w,l,prefix,chrom);
    }
    private static String getPositionsFname(String prefix, int chrom) {
        return prefix + chrom + ".spositions";
    }
    private static String getWeightsFname(String prefix, int chrom) {
        return prefix + chrom + ".sweights";
    }
    private static String getLaSFname(String prefix, int chrom) {
        return prefix + chrom + ".slas";
    }

}