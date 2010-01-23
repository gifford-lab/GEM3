package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * Represents a list of sorted reads on disk
 */
public class PairedHits extends Hits {

    /* 
     * true if the positions (and chrom) for these PairedHits are
     * the left side, meaning that the other chrom and position are for the right side.
     * If false, then the positions are for the right side and the other chrom and
     * position are for the left side.
     */
    private boolean isLeft;
    /* stores the chromosome of the other side of the pair  */
    private IntBP chroms;
    /* stores the position of the other side of the pair */
    private IntBP otherPositions;

    public PairedHits(String prefix, int chrom, boolean isLeft) throws FileNotFoundException, SecurityException, IOException {
        super(chrom,
              getPositionsFname(prefix,chrom,isLeft),
              getWeightsFname(prefix,chrom,isLeft), 
              getLaSFname(prefix,chrom,isLeft));
        this.isLeft = isLeft;
        chroms = openIntBP(getChromsFname(prefix,chrom,isLeft));
        otherPositions = openIntBP(getOtherPosFname(prefix,chrom,isLeft));
    }

    public static void writePairedHits(IntBP positions,
                                       FloatBP weights,
                                       IntBP las,
                                       IntBP otherChroms,
                                       IntBP otherPositions,
                                       String prefix,
                                       int chrom,
                                       boolean isLeft) throws IOException {
        String postmp = getPositionsFname(prefix,chrom,isLeft) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom,isLeft) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom,isLeft) + ".tmp";
        String chrtmp = getChromsFname(prefix,chrom,isLeft) + ".tmp";
        String optmp = getOtherPosFname(prefix,chrom,isLeft) + ".tmp";
        RandomAccessFile positionsRAF = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(lastmp,"rw");
        RandomAccessFile chromsRAF = new RandomAccessFile(chrtmp,"rw");
        RandomAccessFile otherposRAF = new RandomAccessFile(optmp,"rw");

        Bits.sendBytes(positions.bb, 0, positions.bb.limit(), positionsRAF.getChannel());
        Bits.sendBytes(weights.bb, 0, weights.bb.limit(), weightsRAF.getChannel());
        Bits.sendBytes(las.bb, 0, las.bb.limit(), lasRAF.getChannel());
        Bits.sendBytes(otherChroms.bb, 0, otherChroms.bb.limit(), chromsRAF.getChannel());
        Bits.sendBytes(otherPositions.bb, 0, otherPositions.bb.limit(), otherposRAF.getChannel());
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();
        chromsRAF.close();
        otherposRAF.close();

        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(getPositionsFname(prefix,chrom));
        (new File(weightstmp)).renameTo(getWeightsFname(prefix,chrom));
        (new File(lastmp)).renameTo(getLaSFname(prefix,chrom));
        (new File(chrtmp)).renameTo(getChromsFname(prefix,chrom,isLeft));
        (new File(optmp)).renameTo(getOtherPosFname(prefix,chrom,isLeft));        
    }
    public static void writePairedHits(List<PairedHits> hits,
                                       String prefix, 
                                       int chrom,
                                       boolean isLeft) throws IOException {
        Collections.sort(hits, isLeft ? new PairedHitLeftComparator() : new PAiredHitRightComparator());

        IntBP p = new IntBP(hits.length());
        FloatBP w = new FloatBP(hits.length());
        IntBP l = new IntBP(hits.length());
        IntBP c = new IntBP(hits.length());
        Int op = new IntBP(hits.length());
        for (int i = 0; i < hits.length(); i++) {
            PairedHit h = hits.get(i);
            w.put(i, h.weight);
            if (isLeft) {
                p.put(i, h.leftPos);
                l.put(i, makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                c.put(h.rightChrom);
                op.put(h.rightPos);
            } else {
                p.put(i, h.rightPos);
                l.put(i, makeLAS(h.rightLength, h.rightStrand, h.leftLength, h.leftStrand));
                c.put(h.leftChrom);
                op.put(h.leftPos);
            }
        }
        writePairedHits(p,w,l,c,op,prefix,chrom,isLeft);
    }
    private static String getLeftRightSuffix(boolean isLeft) {
        return isLeft ? ".prleft" : ".prright";
    }
    private static String getPositionsFname(String prefix, int chrom, boolean isLeft) {
        prefix + chrom + getLeftRightSuffix(isLeft) + ".positions";
    }
    private static String getWeightsFname(String prefix, int chrom, boolean isLeft) {
        prefix + chrom + getLeftRightSuffix(isLeft) + ".weights";
    }
    private static String getLaSFname(String prefix, int chrom, boolean isLeft) {
        prefix + chrom + getLeftRightSuffix(isLeft) + ".las";
    }
    private static String getChromsFname(String prefix, int chrom, boolean isLeft) {
        prefix + chrom + getLeftRightSuffix(isLeft) + ".chroms";
    }
    private static String getOtherPosFname(String prefix, int chrom, boolean isLeft) {
        prefix + chrom + getLeftRightSuffix(isLeft) + ".otherpos";
    }



}