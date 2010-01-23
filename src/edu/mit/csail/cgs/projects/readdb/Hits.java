package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

public abstract class Hits implements Closeable {

    public static IntBP emptyIntBP = new IntBP(0);
    public static FloatBP emptyFloatBP = new FloatBP(0);
    private static int strandOneMask = 0x00008000;
    private static int lenOneMask = 0x00007FFF;
    private static int strandTwoMask = 0x80000000;
    private static int lenTwoMask = 0x7FFF0000;

    public static IntBP openIntBP(String fname) throws FileNotFoundException, SecurityException, IOException {
        RandomAccessFile raf = null;
        FileChannel fc = null;
        IntBP ib = null;
        ByteBuffer bb = null;
        IOException ioex = null;
        FileNotFoundException fnfex = null;
        SecurityException secex = null;
        try {
            raf = new RandomAccessFile(fname,"r");
            fc = raf.getChannel();
            bb = fc.map(FileChannel.MapMode.READ_ONLY,
                        0,
                        fc.size());
            bb.order(ByteOrder.nativeOrder());
            ib = new IntBuffer(bb);
        } catch (IOException e) {
            ioex = e;
        } catch (FileNotFoundException e) {
            fnfex = e;
        } catch (SecurityException e) {
            secex = e;
        } finally { 
            if (fc != null) {
                fc.close();
            }
            if (raf != null) {
                raf.close();
            }
        }
        if (ioex != null) {
            bb = null;
            ib = null;
            throw ioex;
        }
        if (fnfex != null) {
            bb = null;
            ib = null;
            throw fnfex;
        }
        if (secex != null) {
            bb = null;
            ib = null;
            throw secex;
        }
        return ib;
    }
    public static FloatBP openFloatBP(String fname) throws FileNotFoundException, SecurityException, IOException {
        RandomAccessFile raf = null;
        FileChannel fc = null;
        FloatBP fb = null;
        ByteBuffer bb = null;
        IOException ioex = null;
        FileNotFoundException fnfex = null;
        SecurityException secex = null;
        try {
            raf = new RandomAccessFile(fname,"r");
            fc = raf.getChannel();
            bb = fc.map(FileChannel.MapMode.READ_ONLY,
                        0,
                        fc.size());
            bb.order(ByteOrder.nativeOrder());
            fb = new FloatBuffer(bb);
        } catch (IOException e) {
            ioex = e;
        } catch (FileNotFoundException e) {
            fnfex = e;
        } catch (SecurityException e) {
            secex = e;
        } finally { 
            if (fc != null) {
                fc.close();
            }
            if (raf != null) {
                raf.close();
            }
        }
        if (ioex != null) {
            bb = null;
            fb = null;
            throw ioex;
        }
        if (fnfex != null) {
            bb = null;
            fb = null;
            throw fnfex;
        }
        if (secex != null) {
            bb = null;
            fb = null;
            throw secex;
        }
        return fb;
    }

    /** returns just the length from an int representing a
        length and strand */
    public static short getLengthOne(int las) {
        return (short)(las & lenOneMask);
    }
    public static short getLengthTwo(int las) {
        return (short)((las & lenTwoMask) >> 16);
    }
    public static boolean getStrandOne(int las) {
        return (las & strandOneMask) == 0;
    }
    public static boolean getStrandTwo(int las) {
        return (las & strandTwoMask) == 0;
    }
    public static int makeLAS(short lengthOne, boolean strandOne) {
        return (lengthOne & lenOneMask) | ((strandOne ? 1 : 0) << 15);
    }
    public static int makeLAS(short lengthOne, boolean strandOne, short lengthTwo, boolean strandTwo) {
        return (lengthOne & lenOneMask) | ((strandOne ? 1 : 0) << 15) |
            ((lengthTwo & lenOneMask) << 16) | ((strandTwo ? 1 : 0) << 31);
    }

    private IntBP positions;
    private FloatBP weights;       
    private IntBP lenAndStrand;
    private int chrom;

    public Hits (int chrom, String positionsFname, String weightsFname, String lasFname) throws FileNotFoundException, SecurityException, IOException {
        this.chrom = chrom;
        positions = openIntBP(positionsFname);
        weights = openFloatBP(weightsFname);
        lenAndStrand = openIntBP(lasFname);
    }
    /** gets the buffer of positions */
    public IntBP getPositionsBuffer() {
        return positions;
    }
    /** returns the buffer of weights */
    public FloatBP getWeightsBuffer() {
        return weights;
    }
    /** returns the buffer of lengths and strands */
    public IntBP getLASBuffer() {
        return lenAndStrand
    }
    /**
     * returns indices = int[2] 
     * such that indices[0] is the first element of positions >= startpos
     * and indices[1] is the first element of positions > lastpos.
     * indices[0] == indices[1] -> no hits between startpos and lastpos
     * firstindex is an inclusive lower bound on the index of startpos.
     * lastindex is an inclusive upper bound on the index of lastpos
     */
    public int[] getIndices(int firstindex, int lastindex, int startpos, int lastpos) {
        assert(startpos <= lastpos);
        assert(firstindex <= lastindex);
        assert(firstindex >= 0);
        assert(lastindex >= firstindex);
        assert(lastindex < positions.ib.limit());
        //        System.err.println(String.format("%s Looking for %d - %d in %d - %d", this.toString(), startpos,lastpos,firstindex,lastindex));
        int indices[] = new int[2]; // will be the output
        while (firstindex < positions.limit() && positions.get(firstindex) < startpos) {
            //            System.err.println(String.format("first positions[%d] = %d", firstindex, positions.get(firstindex)));
            firstindex++;
        }
        while (lastindex > firstindex && positions.get(lastindex - 1) > lastpos) {
            //            System.err.println(String.format("last positions[%d] = %d", lastindex, positions.get(lastindex)));
            lastindex--;
        }
        indices[0] = firstindex;
        indices[1] = lastindex;
        //         System.err.println(String.format("Returning %d - %d (%d,%d) = %d", firstindex, lastindex, 
        //                                          positions.get(firstindex), positions.get(lastindex),lastindex - firstindex));
        return indices;

    }
    /**
     * Returns the number of elements between start and stop
     * where firstindex and lastindex are the lower and upper bounds to search
     */
    public int getCountBetween (int firstindex,
                                int lastindex,
                                int start,
                                int stop) throws IOException {       
        int[] p = getIndices(firstindex, lastindex, start,stop);
        return p[1] - p[0];
    }
    public int getCountBetween (int firstindex,
                                int lastindex,
                                int start,
                                int stop,
                                float minweight) throws IOException {       
        int[] p = getIndices(firstindex, lastindex, start,stop);
        int count = 0;
        for (int i = p[0]; i < p[1]; i++) {
            if (weights.get(i) >= minweight) {
                count++;
            }
        }
        return count;
    }
    private IntBP getIntsBetween(IntBP buffer,
                                 int firstindex,
                                 int lastindex,
                                 int start,
                                 int stop) throws IOException {
        int[] p = getIndices(firstindex, lastindex, start,stop);
        if (p[0] >= p[1]) {
            return emptyIntBP;
        }
        return buffer.slice(p[0], p[1] - p[0]);
    }
    /** returns the hits between firstindex and last index and between start and stop, inclusive.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    public IntBP getHitsBetween(int firstindex,
                                int lastindex,
                                int start,
                                int stop) throws IOException {
        return getIntsBetween(positions,firstindex,lastindex,start,stop);
    }
    public IntBP getLASBetween(int firstindex,
                                int lastindex,
                                int start,
                                int stop) throws IOException {
        return getIntsBetween(lenAndStrand,firstindex,lastindex,start,stop);
    }
    /** returns the ints between firstindex and last index and between start and stop, inclusive,
     * for which the weight is greater than minweight.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    private IntBP getIntsBetween(IntBP buffer,
                                 int firstindex,
                                 int lastindex,
                                 int start,
                                 int stop,
                                 float minweight) throws IOException {
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        int n = 0;
        for (int i = p[0]; i < p[1]; i++) {
            if (weights.get(i) >= minweight) {
                n++;
            }
        }
        if (n == 0) {
            return emptyIntBP;
        }

        IntBP output = new IntBP(ByteBuffer.allocate(n*4));
        n = 0;
        for (int i = p[0]; i < p[1]; i++) {
            if (weights.get(i) >= minweight) {
                output.ib.put(n, buffer.get(i));
            }
        }
        return output;
    }
    public IntBP getHitsBetween(IntBP buffer,
                                 int firstindex,
                                 int lastindex,
                                 int start,
                                 int stop,
                                 float minweight) throws IOException {
        return getIntsBetween(positions,firstindex,lastindex,start,stop,minweight);
    }
    public IntBP getLASBetween(IntBP buffer,
                               int firstindex,
                               int lastindex,
                               int start,
                               int stop,
                               float minweight) throws IOException {
        return getIntsBetween(lenAndStrand,firstindex,lastindex,start,stop,minweight);
    }
    /** returns the weights between firstindex and last index and between start and stop, inclusive.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    public FloatBP getWeightsBetween(int firstindex,
                                     int lastindex,
                                     int start,
                                     int stop) throws IOException {
        int[] p = getIndices(firstindex, lastindex, start,stop);
        if (p[0] >= p[1]) {
            return emptyFloatBP;
        }
        return weights.slice(p[0], p[1] - p[0]);
    }
    /** returns the weights between firstindex and last index and between start and stop, inclusive,
     * for which the weight is greater than minweight.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    public FloatBP getWeightsBetween(int firstindex,
                                     int lastindex,
                                     int start,
                                     int stop,
                                     float minweight) throws IOException {
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        int n = 0;
        for (int i = p[0]; i < p[1]; i++) {
            if (weights.get(i) >= minweight) {
                n++;
            }
        }
        if (n == 0) {
            return emptyFloatBP;
        }
        FloatBP output = new FloatBP(ByteBuffer.allocate(n*4));
        n = 0;
        for (int i = p[0]; i < p[1]; i++) {
            if (weights.get(i) >= minweight) {
                output.fb.put(n,weights.get(i));
            }
        }
        return output;
    }
    /** return a histogram from positions start to stop
     *  in units of stepsize.  firstindex and lastindex come from
     *  Header.getFirstIndex(start) and Header.getLastIndex(stop)
     *
     *  If extension is non-zero, then each hit is extended to
     * cover a region of extension bp (positive values extend to greater
     * coordinates, negative numbers to smaller coordinates) and the
     * hit is counted in each bin that it touches.
     */                                           
    public int[] histogram(int firstindex,
                           int lastindex,
                           int start,
                           int stop,
                           int stepsize,
                           boolean extension) throws IOException {
        IntBP positions = getPositionsBuffer();
        int output[] = new int[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (!extension) {
            for (int i = p[0]; i < p[1]; i++) {
                output[(positions.get(i) - start) / stepsize]++;            
            }
        } else {
            IntBP las = getLASBuffer();
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                int l = las.get(i);
                short len = getLengthOne(l);
                boolean strand = getStrandOne(l);
                if (strand) {
                    while (len > 0 && ++bin < output.length) {
                        output[bin]++;
                        len -= stepsize;
                    }
                } else {
                    while (len > 0 && --bin > 0) {
                        output[bin]++;
                        len -= stepsize;
                    }
                }
            }
        }
        return output;
    }
    /**
     * histogram with a minimum weight
     */
    public int[] histogram(int firstindex,
                           int lastindex,
                           int start,
                           int stop,
                           int stepsize,
                           float minweight,
                           boolean extension) throws IOException {
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
        int output[] = new int[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (!extension) {
            for (int i = p[0]; i < p[1]; i++) {
                if (weights.get(i) >= minweight) {
                    output[(positions.get(i) - start) / stepsize]++;            
                }
            }
        } else {
            IntBP las = getLASBuffer();
            for (int i = p[0]; i < p[1]; i++) {
                if (weights.get(i) >= minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    int l = las.get(i);
                    short len = getLengthOne(l);
                    boolean strand = getStrandOne(l);
                    if (strand) {
                        while (len > 0 && ++bin < output.length) {
                            output[bin]++;
                            len -= stepsize;
                        }
                    } else {
                        while (len > 0 && --bin > 0) {
                            output[bin]++;
                            len -= stepsize;
                        }
                    }
                }
            }
        }
        return output;
    }
    public float[] weightHistogram(int firstindex,
                                   int lastindex,
                                   int start,
                                   int stop,
                                   int stepsize,
                                   boolean extension) throws IOException {
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
        float output[] = new int[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (!extension) {
            for (int i = p[0]; i < p[1]; i++) {
                output[(positions.get(i) - start) / stepsize] += weights.get(i);            
            }
        } else {
            IntBP las = getLASBuffer();
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                int l = las.get(i);
                short len = getLengthOne(l);
                boolean strand = getStrandOne(l);
                if (strand) {
                    while (len > 0 && ++bin < output.length) {
                        output[bin] += weights.get(i);
                        len -= stepsize;
                    }
                } else {
                    while (len > 0 && --bin > 0) {
                        output[bin] += weights.get(i);
                        len -= stepsize;
                    }
                }
            }
        }
        return output;
    }
    /**
     * histogram with a minimum weight
     */
    public float[] weightHistogram(int firstindex,
                                   int lastindex,
                                   int start,
                                   int stop,
                                   int stepsize,
                                   float minweight,
                                   boolean extension) throws IOException {
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
        float output[] = new float[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (!extension) {
            for (int i = p[0]; i < p[1]; i++) {
                float f = weights.get(i);
                if (f > minweight) {
                    output[(positions.get(i) - start) / stepsize] += f;
                }
            }
        } else {
            IntBP las = getLASBuffer();
            for (int i = p[0]; i < p[1]; i++) {
                float f = weights.get(i);
                if (f >= minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    int l = las.get(i);
                    short len = getLengthOne(l);
                    boolean strand = getStrandOne(l);
                    if (strand) {
                        while (len > 0 && ++bin < output.length) {
                            output[bin] += f;
                            len -= stepsize;
                        }
                    } else {
                        while (len > 0 && --bin > 0) {
                            output[bin] += f;
                            len -= stepsize;
                        }
                    }
                }
            }
        }
        return output;
    }   
    
    public void close() throws IOException {        
        positions.ib = null;
        positions.bb = null;
        positions = null;
        weights.fb = null;
        weights.bb = null;
        weights = null;
        lenAndStrand.ib = null;
        lenAndStrand.bb = null;
        lenAndStrand = null;
    }


}