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
    /** appends a sorted list of hits to an existing set of hits. */
    public void appendSingleHits(SingleHit[] hits,
                                 String prefix,
                                 int chrom) throws IOException {
        if (hits.length == 0) {
            return;
        }
        int s = getPositionsBuffer().limit() - 1;
        SingleHit last = new SingleHit(chrom, getPositionsBuffer().get(s), getWeightsBuffer().get(s),
                                       Hits.getStrandOne(getLASBuffer().get(s)),
                                       Hits.getLengthOne(getLASBuffer().get(s)));
        if (hits[0].compareTo(last) > 0 ) {
            append(hits,prefix,chrom);
        } else {
            merge(hits,prefix,chrom);
        }
    }
    private void append(SingleHit[] hits,
                       String prefix,
                       int chrom) throws IOException {
        RandomAccessFile  positionsRAF = new RandomAccessFile(getPositionsFname(prefix,chrom),"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(getWeightsFname(prefix,chrom),"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(getLaSFname(prefix,chrom),"rw");
        positionsRAF.seek(positionsRAF.length());
        weightsRAF.seek(weightsRAF.length());
        lasRAF.seek(lasRAF.length());
        for (int i = 0; i < hits.length; i++) {
            SingleHit h = hits[i];
            positionsRAF.writeInt(h.pos);
            weightsRAF.writeFloat(h.weight);
            lasRAF.writeInt(makeLAS(h.length, h.strand));
        }
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();        
    }
    private void merge(SingleHit[] hits,
                       String prefix,
                       int chrom) throws IOException {
        String postmp = getPositionsFname(prefix,chrom) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom) + ".tmp";
        RandomAccessFile positionsRAF = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(lastmp,"rw");
        int newsize = getPositionsBuffer().limit() + hits.length;
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
        IntBP oldpositions = getPositionsBuffer();
        FloatBP oldweights = getWeightsBuffer();
        IntBP oldlas = getLASBuffer();
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
        posfile = null;
        weightfile = null;
        lasfile = null;        
        oldpositions = null;
        oldweights =null;
        oldlas = null;
        positionsRAF.close();
        weightsRAF.close();
        lasRAF.close();
        /* ideally this part with the renames would atomic... */
        (new File(postmp)).renameTo(new File(getPositionsFname(prefix,chrom)));
        (new File(weightstmp)).renameTo(new File(getWeightsFname(prefix,chrom)));
        (new File(lastmp)).renameTo(new File(getLaSFname(prefix,chrom)));
    }
    public void resort(String prefix, int chrom) throws IOException {
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
        IntBP las = getLASBuffer();
        long indices[] = new long[positions.limit()];
        for (int i = 0; i < indices.length; i++) {
            long v = positions.get(i);
            v <<= 32;
            v |= i;
            indices[i] = v;
        }
        Arrays.sort(indices);

        String postmp = getPositionsFname(prefix,chrom) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom) + ".tmp";
        RandomAccessFile positionsRAF = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weightsRAF = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile lasRAF = new RandomAccessFile(lastmp,"rw");
        int newsize = getPositionsBuffer().limit();
        IntBP posfile = new IntBP(positionsRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                                0,
                                                                newsize * 4));
        FloatBP weightfile = new FloatBP(weightsRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                                     0,
                                                                     newsize * 4));
        IntBP lasfile = new IntBP(lasRAF.getChannel().map(FileChannel.MapMode.READ_WRITE,
                                                          0,
                                                          newsize * 4));

        for (int i = 0; i < indices.length; i++) {
            int index = (int)(indices[i] & 0xffffffffL);
            int pos = (int)(indices[i] >> 32);
            posfile.put(i, pos);
            weightfile.put(i, weights.get(index));
            lasfile.put(i, las.get(index));
        }
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