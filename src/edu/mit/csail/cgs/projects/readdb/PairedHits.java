package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.List;
import java.util.Collections;
import java.util.Comparator;
import java.util.Arrays;

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
    public boolean isLeft () {return isLeft;}
    public IntBP getChromsBuffer() {return chroms;}
    public IntBP getOtherPosBuffer() {return otherPositions;}
    public IntBP getOtherChromsBetween(int firstindex,
                                       int lastindex,
                                       int start,
                                       int stop,
                                       Float minweight,
                                       Boolean isPlus) throws IOException {
        return getIntsBetween(chroms,firstindex,lastindex,start,stop,minweight,isPlus);
    }
    public IntBP getOtherPositionsBetween(int firstindex,
                                          int lastindex,
                                          int start,
                                          int stop,
                                          Float minweight,
                                          Boolean isPlus) throws IOException {
        return getIntsBetween(otherPositions,firstindex,lastindex,start,stop,minweight,isPlus);
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
        if (otherChroms.bb.order() != ByteOrder.nativeOrder()) {
            Bits.flipByteOrder(otherChroms.ib);
            otherChroms.bb.order(ByteOrder.nativeOrder());
        }
        if (otherPositions.bb.order() != ByteOrder.nativeOrder()) {
            Bits.flipByteOrder(otherPositions.ib);
            otherPositions.bb.order(ByteOrder.nativeOrder());
        }

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
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom,isLeft)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom,isLeft)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom,isLeft)));
        (new File(chrtmp)).renameTo(new File(getChromsFname(prefix,chrom,isLeft)));
        (new File(optmp)).renameTo(new File(getOtherPosFname(prefix,chrom,isLeft)));        
    }
    public static void writePairedHits(PairedHit[] hits,
                                       String prefix, 
                                       int chrom,
                                       boolean isLeft) throws IOException {
        Comparator<PairedHit> comp = isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator();
        boolean sorted = true;
        int i = 1;
        while (sorted && i < hits.length) {
            if (comp.compare(hits[i-1],hits[i]) >= 0) {
                sorted = false;
            }
            i++;
        }
        if (!sorted) {
            Arrays.sort(hits, comp);
        }
        //        System.err.println("STORING HITS " + hits);
        IntBP p = new IntBP(hits.length);
        FloatBP w = new FloatBP(hits.length);
        IntBP l = new IntBP(hits.length);
        IntBP c = new IntBP(hits.length);
        IntBP op = new IntBP(hits.length);
        for (i = 0; i < hits.length; i++) {
            PairedHit h = hits[i];
            w.put(i, h.weight);
            if (isLeft) {
                p.put(i, h.leftPos);
                l.put(i, makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                c.put(i, h.rightChrom);
                op.put(i, h.rightPos);
            } else {
                p.put(i, h.rightPos);
                l.put(i, makeLAS(h.rightLength, h.rightStrand, h.leftLength, h.leftStrand));
                c.put(i, h.leftChrom);
                op.put(i, h.leftPos);
            }
        }
        writePairedHits(p,w,l,c,op,prefix,chrom,isLeft);
    }
    private static String getLeftRightSuffix(boolean isLeft) {
        return isLeft ? ".prleft" : ".prright";
    }
    private static String getPositionsFname(String prefix, int chrom, boolean isLeft) {
        return prefix + chrom + getLeftRightSuffix(isLeft) + ".positions";
    }
    private static String getWeightsFname(String prefix, int chrom, boolean isLeft) {
        return prefix + chrom + getLeftRightSuffix(isLeft) + ".weights";
    }
    private static String getLaSFname(String prefix, int chrom, boolean isLeft) {
        return prefix + chrom + getLeftRightSuffix(isLeft) + ".las";
    }
    private static String getChromsFname(String prefix, int chrom, boolean isLeft) {
        return prefix + chrom + getLeftRightSuffix(isLeft) + ".chroms";
    }
    private static String getOtherPosFname(String prefix, int chrom, boolean isLeft) {
        return prefix + chrom + getLeftRightSuffix(isLeft) + ".otherpos";
    }



}