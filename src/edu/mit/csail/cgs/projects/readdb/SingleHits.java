package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.List;
import java.util.Collections;
import java.util.Arrays;

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
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom)));
    }
    public static void writeSingleHits(SingleHit[] hits,
                                       String prefix, 
                                       int chrom) throws IOException {
        IntBP p = new IntBP(hits.length);
        FloatBP w = new FloatBP(hits.length);
        IntBP l = new IntBP(hits.length);
        for (int i = 0; i < hits.length; i++) {
            SingleHit h = hits[i];
            p.put(i, h.pos);
            w.put(i, h.weight);
            l.put(i, makeLAS(h.length, h.strand));
        }
        writeSingleHits(p,w,l,prefix,chrom);
    }
    public static void appendSingleHits(SingleHit[] hits,
                                        SingleHits oldhits,
                                        String prefix,
                                        int chrom) throws IOException {
        IntBP oldpositions = oldhits.getPositionsBuffer();
        FloatBP oldweights = oldhits.getWeightsBuffer();
        IntBP oldlas = oldhits.getLASBuffer();

        int newsize = oldpositions.limit() + hits.length;
        
        String postmp = getPositionsFname(prefix,chrom) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom) + ".tmp";
        RandomAccessFile positionsRAF = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(lastmp,"rw");
        IntBP posfile = new IntBP(positionsRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                                0,
                                                                newsize * 4));
        FloatBP weightfile = new FloatBP(weightsRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                                     0,
                                                                     newsize * 4));
        IntBP lasfile = new IntBP(lasRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                          0,
                                                          newsize * 4));
        
        int oldp = 0;
        int newp = 0;
        int pos = 0;
        //        System.err.println(String.format("APPENDING for %d : %d and %d",chrom, oldpositions.limit(), hits.length));
        while (oldp < oldpositions.limit() || newp < hits.length) {
            while (newp < hits.length && (oldp == oldpositions.limit() || hits[newp].pos <= oldpositions.get(oldp))) {
                posfile.put(pos, hits[newp].pos);
                weightfile.put(pos, hits[newp].weight);
                lasfile.put(pos, Hits.makeLAS(hits[newp].length, hits[newp].strand));
                newp++;
                pos++;
            }
            while (oldp < oldpositions.limit() && (newp == hits.length || oldpositions.get(oldp) <= hits[newp].pos)) {
                posfile.put(pos,oldpositions.get(oldp));
                weightfile.put(pos,oldweights.get(oldp));
                lasfile.put(pos,oldlas.get(oldp));
                oldp++;
                pos++;                
            }          
            //            System.err.println(String.format("%d %d %d", pos, newp, oldp));
        }
        oldpositions = null;
        oldweights =null;
        oldlas = null;
        posfile = null;
        weightfile = null;
        lasfile = null;
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();
        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom)));
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