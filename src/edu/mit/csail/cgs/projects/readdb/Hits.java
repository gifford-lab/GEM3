package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/** 
 * Represents the list of sorted read positions on disk
 */
public class Hits implements Closeable {

    private IntBP positions;
    private FloatBP weights;       
    private static IntBP emptyIntBP = new IntBP(0);
    private static FloatBP emptyFloatBP = new FloatBP(0);

    /**
     * Initializes a Hits object from a file
     */
    public Hits (String positionsfname, String weightsfname) throws FileNotFoundException, SecurityException, IOException {
        try {
            RandomAccessFile positionsRAF = null, weightsRAF = null;
            FileChannel positionsFC = null, weightsFC = null;
            IOException ioex = null;
            ByteBuffer bb;
            try {
                positionsRAF = new RandomAccessFile(positionsfname,"r");
                positionsFC = positionsRAF.getChannel();
                bb = positionsFC.map(FileChannel.MapMode.READ_ONLY,
                                     0,
                                     positionsFC.size());
                bb.order(ByteOrder.nativeOrder());
                positions = new IntBP(bb);
                weightsRAF = new RandomAccessFile(weightsfname,"r");
                weightsFC = weightsRAF.getChannel();
                bb = weightsFC.map(FileChannel.MapMode.READ_ONLY,
                                   0,
                                   weightsFC.size());
                bb.order(ByteOrder.nativeOrder());
                weights = new FloatBP(bb);
            } catch (IOException e) {
                ioex = e;
                bb = null;
                positions = null;
                weights = null;
            } finally {
                if (positionsFC != null) { positionsFC.close();  }
                if (positionsRAF != null) { positionsRAF.close();}
                if (weightsFC != null) {weightsFC.close(); }
                if (weightsRAF != null) {weightsRAF.close(); }
            }            
            if (ioex != null) {
                throw ioex;
            }
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        }
    }
    /** 
     * gets the buffer of positions
     */
    public IntBP getPositionsBuffer() {
        return positions;
    }
    /** returns the buffer of weights
     */
    public FloatBP getWeightsBuffer() {return weights;}
    /**
     * returns indices = int[2] 
     * such that indices[0] is the first element of positions >= startpos
     * and indices[1] is the first element of positions > lastpos.
     * indices[0] == indices[1] -> no hits between startpos and lastpos
     * firstindex is an inclusive lower bound on the index of startpos.
     * lastindex is an inclusive upper bound on the index of lastpos
     */
    public int[] getIndicesLinear(int firstindex, int lastindex, int startpos, int lastpos) {
        assert(startpos <= lastpos);
        assert(firstindex <= lastindex);
        assert(firstindex >= 0);
        assert(lastindex >= firstindex);
        assert(lastindex < positions.ib.limit());
        //        System.err.println(String.format("%s Looking for %d - %d in %d - %d", this.toString(), startpos,lastpos,firstindex,lastindex));
        int indices[] = new int[2]; // will be the output
        int limit = positions.ib.limit();
        if (lastindex < limit) {
            limit = lastindex;
        }

        while (firstindex < limit && positions.ib.get(firstindex) < startpos) {
            //            System.err.println(String.format("first positions[%d] = %d", firstindex, positions.get(firstindex)));
            firstindex++;
        }
        while (lastindex > firstindex && positions.ib.get(lastindex - 1) > lastpos) {
            //            System.err.println(String.format("last positions[%d] = %d", lastindex, positions.get(lastindex)));
            lastindex--;
        }
        indices[0] = firstindex;
        indices[1] = lastindex;
        //         System.err.println(String.format("Returning %d - %d (%d,%d) = %d", firstindex, lastindex, 
        //                                          positions.get(firstindex), positions.get(lastindex),lastindex - firstindex));
        return indices;

    }
    /* don't actually use this version.  You might think that binary search is faster, but it turns
       out not to be.  The linear search is over a limited number of elements, is very simple code,
       and has a great data cache hit rate
    */      
    public int[] getIndices(int firstindex, int lastindex, int startpos, int lastpos) {
        assert(startpos <= lastpos);
        assert(firstindex <= lastindex);
        int indices[] = new int[2]; // will be the output
        /* this is a binary search in two steps.  In the first step, we establish the upper bound for the binary
           search by looking out in windows of doubling size until we've bounded the spot we're looking for.  
           We do this because a binary search over the entire array of hits would defeat the purpose of the
           index in Header.  The second step is the actual binary search using firstindex as the lower bound and the
           upper bound that we found
        */
        int step = 256;
        int limit = positions.ib.limit();
        if (lastindex < limit) {
            limit = lastindex;
        }
        while (firstindex + step < limit &&
               positions.ib.get(firstindex + step) < startpos) {
            firstindex += step;
            step *= 2;
        }
        while (firstindex < limit && 
               positions.ib.get(firstindex) < startpos) {
            step = (step >> 2) | 1;
            if (step == 1 ||
                (firstindex + step < limit &&
                 positions.ib.get(firstindex + step - 1) < startpos)) {
                firstindex += step;
            }
        }

        step = 256;
        while (lastindex - step > firstindex &&
               positions.ib.get(lastindex - step) > lastpos) {
            lastindex -= step;
            step *= 2;
        }
        while (lastindex > firstindex && 
               positions.ib.get(lastindex - 1) > lastpos) {
            step = (step >> 2) | 1;
            if ((step == 1 || lastindex - step > firstindex) &&
                positions.ib.get(lastindex - step) > lastpos) {
                lastindex -= step;
            }
        }

        indices[0] = firstindex;
        indices[1] = lastindex;

        assert(indices[0] <= indices[1]);
        assert(indices[0] >= 0);
        assert(indices[1] <= positions.size());
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
    
    /** returns the hits between firstindex and last index and between start and stop, inclusive.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    public IntBP getHitsBetween(int firstindex,
                                int lastindex,
                                int start,
                                int stop) throws IOException {
        int[] p = getIndices(firstindex, lastindex, start,stop);
        if (p[0] >= p[1]) {
            return emptyIntBP;
        }
        return positions.slice(p[0], p[1] - p[0]);
    }
    /** returns the hits between firstindex and last index and between start and stop, inclusive,
     * for which the weight is greater than minweight.
     *  firstindex and lastindex come from Header.getFirstIndex and Header.getLastIndex
     */
    public IntBP getHitsBetween(int firstindex,
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
                output.ib.put(n, positions.get(i));
            }
        }
        return output;
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
                           int extension) throws IOException {
        extension /= stepsize;
        int output[] = new int[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (extension == 0) {
            for (int i = p[0]; i < p[1]; i++) {
                output[(positions.get(i) - start) / stepsize]++;            
            }
        } else if (extension > 0) {
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                for (int j = 0; j <= extension && bin+j < output.length; j++) {
                    output[bin+j]++;            
                }
            }
        } else {
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                for (int j = 0; j >= extension && bin+j >= 0; j--) {
                    output[bin+j]++;            
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
                           int extension) throws IOException {
        extension /= stepsize;
        int output[] = new int[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (extension == 0) {
            for (int i = p[0]; i < p[1]; i++) {
                if (weights.get(i) >= minweight) {
                    output[(positions.get(i) - start) / stepsize]++;            
                }
            }
        } else if (extension > 0) {
            for (int i = p[0]; i < p[1]; i++) {
                if (weights.get(i) >= minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    for (int j = 0; j <= extension && bin+j < output.length; j++) {
                        output[bin+j]++;
                    }
                }
            }
        } else {
            for (int i = p[0]; i < p[1]; i++) {
                if (weights.get(i) >= minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    for (int j = 0; j >= extension && bin+j >= 0; j--) {
                        output[bin+j]++;
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
                                   int extension) throws IOException {
        extension /= stepsize;
        float output[] = new float[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (extension == 0) {
            for (int i = p[0]; i < p[1]; i++) {
                output[(positions.get(i) - start) / stepsize] += weights.get(i);            
            }
        } else if (extension > 0) {
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                float f = weights.get(i);
                for (int j = 0; j <= extension && bin+j < output.length; j++) {
                    output[bin+j] += f;
                }
            }
        } else {
            for (int i = p[0]; i < p[1]; i++) {
                int bin = (positions.get(i) - start) / stepsize;
                float f = weights.get(i);
                for (int j = 0; j >= extension && bin+j>=0; j--) {
                    output[bin+j] += f;
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
                                   int extension) throws IOException {
        float output[] = new float[(stop - start) / stepsize + 1];
        for (int j = 0; j < output.length; j++) {
            output[j] = 0;
        }
        int[] p = getIndices(firstindex, lastindex, start,stop);        
        if (extension == 0) {
            for (int i = p[0]; i < p[1]; i++) {
                float f = weights.get(i);
                if (f > minweight) {
                    output[(positions.get(i) - start) / stepsize] += f;
                }
            }
        } else if (extension > 0) {
            for (int i = p[0]; i < p[1]; i++) {
                float f = weights.get(i);
                if (f > minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    for (int j = 0; j <= extension && bin+j < output.length; j++) {
                        output[bin+j] += f;
                    }
                }
            }
        } else {
            for (int i = p[0]; i < p[1]; i++) {
                float f = weights.get(i);
                if (f > minweight) {
                    int bin = (positions.get(i) - start) / stepsize;
                    for (int j = 0; j >= extension && bin+j >= 0; j--) {
                        output[bin+j] += f;            
                    }
                }
            }
        }
        return output;
    }   
    
    public void close() throws IOException {        
        positions.ib = null;
        positions.bb = null;
        weights.fb = null;
        weights.bb = null;
        positions = null;
        weights = null;
    }

}