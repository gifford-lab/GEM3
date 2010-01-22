package edu.mit.csail.cgs.projects.readdb;

import java.io.*;

/** 
 * Represents the list of sorted reads on disk
 */
public class SingleHits extends Hits {
    /**
     * Initializes a Hits object from a file
     */
    public SingleHits (String prefix, int chrom) throws FileNotFoundException, SecurityException, IOException {
        try {
            RandomAccessFile positionsRAF = null, weightsRAF = null, lasRAF = null;
            FileChannel positionsFC = null, weightsFC = null, lasFS = null;

            this.chrom = chrom;
            

            try {
                String fname = getPositionsFname(prefix,chrom);

            } catch (IOException e) {
                ioex = e;
            } finally {
                if (positionsFC != null) { positionsFC.close();  }
                if (positionsRAF != null) { positionsRAF.close();}
                if (weightsFC != null) {weightsFC.close(); }
                if (weightsRAF != null) {weightsRAF.close(); }
                if (lasRAF != null) {lasRAF.close();}
                if (lasFC != null) {lasFC.close();}
            }            
            if (ioex != null) {
                throw ioex;
            }
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        }
    }
    public static void writeSingleHits(List<SingleHits> hits,
                                       String prefix, 
                                       int chrom) throws IOException {
        Collections.sort(hits);
        String postmp = getPositionsFname(prefix,chrom) + ".tmp";
        String weightstmp = getWeightsFname(prefix,chrom) + ".tmp";
        String lastmp = getLaSFname(prefix,chrom) + ".tmp";

        IntBP p = new IntBP(hits.length());
        FloatBP w = new FloatBP(hits.length());
        IntBP l = new IntBP(hits.length());
        for (int i = 0; i < hits.length(); i++) {
            SingleHit h = hits.get(i);
            p.put(i, h.pos);
            w.put(i, h.weight);
            l.put(i, makeLAS(h.length, h.strand));
        }
        RandomAccessFile positions = new RandomAccessFile(postmp,"rw");
        RandomAccessFile weights = new RandomAccessFile(weightstmp,"rw");
        RandomAccessFile las = new RandomAccessFile(lastmp,"rw");

        Bits.sendBytes(p.bb, 0, hits.length(), positions.getChannel());
        Bits.sendBytes(w.bb, 0, hits.length(), weights.getChannel());
        Bits.sendBytes(l.bb, 0, hits.length(), lastmp.getChannel());
        positions.close();
        weights.close();
        las.close();
    }
    private static String getPositionsFname(String prefix, int chrom) {
        prefix + "." + chrom + ".spositions";
    }
    private static String getWeightsFname(String prefix, int chrom) {
        prefix + "." + chrom + ".sweights";
    }
    private static String getLaSFname(String prefix, int chrom) {
        prefix + "." + chrom + ".slas";
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
        IntBP positions = getPositionsBuffer();
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
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
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
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
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
        IntBP positions = getPositionsBuffer();
        FloatBP weights = getWeightsBuffer();
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
    

}