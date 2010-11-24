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
    private static DataOutputStream dos(String f) throws IOException {
        return new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f)));
    }
    /** hits is a sorted list.
        prefix is file name prefix 
        chrom is chromosome number (must match left if isLeft is true and right otherwise)
    */
    public static void writePairedHits(List<PairedHit> hits,
                                       String prefix, 
                                       int chrom,
                                       boolean isLeft) throws IOException {
        String postmp = getPositionsFname(prefix,chrom,isLeft) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom,isLeft) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom,isLeft) + ".tmp";
        String chrtmp = getChromsFname(prefix,chrom,isLeft) + ".tmp";
        String optmp = getOtherPosFname(prefix,chrom,isLeft) + ".tmp";
        DataOutputStream positionsf = dos(postmp);
        DataOutputStream weightsf = dos(weightstmp);
        DataOutputStream lasf = dos(lastmp);
        DataOutputStream chrf = dos(chrtmp);
        DataOutputStream opf = dos(optmp);

        if (isLeft) {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                positionsf.writeInt(h.leftPos);
                weightsf.writeFloat(h.weight);
                lasf.writeInt(makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                chrf.writeInt(h.rightChrom);
                opf.writeInt(h.rightPos);
            }
        } else {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                positionsf.writeInt(h.rightPos);
                weightsf.writeFloat(h.weight);
                lasf.writeInt(makeLAS(h.rightLength, h.rightStrand, h.leftLength, h.leftStrand));
                chrf.writeInt(h.leftChrom);
                opf.writeInt(h.leftPos);
            }
        }
        positionsf.close();
        weightsf.close();
        lasf.close();
        chrf.close();
        opf.close();
        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom,isLeft)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom,isLeft)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom,isLeft)));
        (new File(chrtmp)).renameTo(new File(getChromsFname(prefix,chrom,isLeft)));
        (new File(optmp)).renameTo(new File(getOtherPosFname(prefix,chrom,isLeft)));                
    }
    public void appendPairedHits(List<PairedHit> hits,
                                 String prefix, 
                                 int chrom,
                                 boolean isLeft) throws IOException {
        if (hits.size() == 0) {
            return;
        }
        int newsize = getPositionsBuffer().limit() + hits.size();
        int s = getPositionsBuffer().limit() - 1;
        PairedHit last = new PairedHit(isLeft ? chrom : chroms.get(s),
                                       isLeft ? getPositionsBuffer().get(s) : otherPositions.get(s),
                                       isLeft ? getStrandOne(getLASBuffer().get(s)) : getStrandTwo(getLASBuffer().get(s)),
                                       isLeft ? getLengthOne(getLASBuffer().get(s)) : getLengthTwo(getLASBuffer().get(s)),

                                       !isLeft ? chrom : chroms.get(s),
                                       !isLeft ? getPositionsBuffer().get(s) : otherPositions.get(s),
                                       !isLeft ? getStrandOne(getLASBuffer().get(s)) : getStrandTwo(getLASBuffer().get(s)),
                                       !isLeft ? getLengthOne(getLASBuffer().get(s)) : getLengthTwo(getLASBuffer().get(s)),
                                       
                                       getWeightsBuffer().get(s));
        Comparator<PairedHit> comparator = isLeft? new PairedHitLeftComparator() : new PairedHitRightComparator();
        if (comparator.compare(hits.get(0), last) > 0) {
            append(hits,prefix,chrom,isLeft);
        } else {
            merge(hits,prefix,chrom,isLeft);
        }
    }
    private void append(List<PairedHit> hits,
                                 String prefix, 
                                 int chrom,
                                 boolean isLeft) throws IOException {
        RandomAccessFile positionsRAF = new RandomAccessFile(getPositionsFname(prefix,chrom,isLeft),"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(getWeightsFname(prefix,chrom,isLeft),"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(getLaSFname(prefix,chrom,isLeft),"rw");
        RandomAccessFile chromsRAF = new RandomAccessFile(getChromsFname(prefix,chrom,isLeft),"rw");
        RandomAccessFile otherposRAF = new RandomAccessFile(getOtherPosFname(prefix,chrom,isLeft),"rw");
        positionsRAF.seek(positionsRAF.length());
        weightsRAF.seek(weightsRAF.length());
        lasRAF.seek(lasRAF.length());
        chromsRAF.seek(chromsRAF.length());
        otherposRAF.seek(otherposRAF.length());
        if (isLeft) {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                positionsRAF.writeInt(h.leftPos);
                weightsRAF.writeFloat(h.weight);
                lasRAF.writeInt(makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                chromsRAF.writeInt(h.rightChrom);
                otherposRAF.writeInt(h.rightPos);

            }
        } else {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                positionsRAF.writeInt(h.leftPos);
                weightsRAF.writeFloat(h.weight);
                lasRAF.writeInt(makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                chromsRAF.writeInt(h.rightChrom);
                otherposRAF.writeInt(h.rightPos);
            }
        }
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();
        chromsRAF.close();
        otherposRAF.close();
    }
    private void merge(List<PairedHit> hits,
                        String prefix, 
                        int chrom,
                        boolean isLeft) throws IOException {
        String postmp = getPositionsFname(prefix,chrom,isLeft) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom,isLeft) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom,isLeft) + ".tmp";
        String chrtmp = getChromsFname(prefix,chrom,isLeft) + ".tmp";
        String optmp = getOtherPosFname(prefix,chrom,isLeft) + ".tmp";
        DataOutputStream positionsf = dos(postmp);
        DataOutputStream weightsf = dos(weightstmp);
        DataOutputStream lasf = dos(lastmp);
        DataOutputStream chrf = dos(chrtmp);
        DataOutputStream opf = dos(optmp);

        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
        IntBP lenAndStrand = getLASBuffer();

        int oldpos = 0;
        if (isLeft) {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                while (oldpos < positions.limit() && 
                       (positions.get(oldpos) < h.leftPos ||
                        (positions.get(oldpos) == h.leftPos &&
                         Hits.getLengthOne(lenAndStrand.get(oldpos)) < h.leftLength))) {
                    positionsf.writeInt(positions.get(oldpos));
                    weightsf.writeFloat(weights.get(oldpos));
                    lasf.writeInt(lenAndStrand.get(oldpos));
                    chrf.writeInt(chroms.get(oldpos));
                    opf.writeInt(otherPositions.get(oldpos));
                    oldpos++;
                }
                positionsf.writeInt(h.leftPos);
                weightsf.writeFloat(h.weight);
                lasf.writeInt(makeLAS(h.leftLength, h.leftStrand, h.rightLength, h.rightStrand));
                chrf.writeInt(h.rightChrom);
                opf.writeInt(h.rightPos);
            }
        } else {
            for (int i = 0; i < hits.size(); i++) {
                PairedHit h = hits.get(i);
                while (oldpos < positions.limit() && 
                       (positions.get(oldpos) < h.rightPos ||
                        (positions.get(oldpos) == h.rightPos &&
                         Hits.getLengthOne(lenAndStrand.get(oldpos)) < h.rightLength))) {
                    positionsf.writeInt(positions.get(oldpos));
                    weightsf.writeFloat(weights.get(oldpos));
                    lasf.writeInt(lenAndStrand.get(oldpos));
                    chrf.writeInt(chroms.get(oldpos));
                    opf.writeInt(otherPositions.get(oldpos));
                    oldpos++;
                }                
                positionsf.writeInt(h.rightPos);
                weightsf.writeFloat(h.weight);
                lasf.writeInt(makeLAS(h.rightLength, h.rightStrand, h.leftLength, h.leftStrand));
                chrf.writeInt(h.leftChrom);
                opf.writeInt(h.leftPos);
            }
        }
        while (oldpos < positions.limit()) {
            positionsf.writeInt(positions.get(oldpos));
            weightsf.writeFloat(weights.get(oldpos));
            lasf.writeInt(lenAndStrand.get(oldpos));
            chrf.writeInt(chroms.get(oldpos));
            opf.writeInt(otherPositions.get(oldpos));
            oldpos++;
        }

        positionsf.close();
        weightsf.close();
        lasf.close();
        chrf.close();
        opf.close();
        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom,isLeft)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom,isLeft)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom,isLeft)));
        (new File(chrtmp)).renameTo(new File(getChromsFname(prefix,chrom,isLeft)));
        (new File(optmp)).renameTo(new File(getOtherPosFname(prefix,chrom,isLeft)));                               
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