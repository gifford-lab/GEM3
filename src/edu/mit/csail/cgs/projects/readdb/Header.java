package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.util.*;
import java.nio.IntBuffer;
import java.nio.ByteOrder;

/**
 * index information for a set of hits
 */
public class Header implements Closeable {

    private int numHits;    
    /* indexPositions and indexPointers are paired arrays.  
       indexPositions contains values from the hits file;
       the corresponding entry in indexPointers is the offset (in sizeof(int) units)
       of the first position in the hits file that contains that value.  
       indexPositions may not contain duplicate values.
    */
    private int[] indexPositions, indexPointers;

    public int getNumHits() {return numHits;}
    /* returns the first element in the file of hits
       that should be read when looking for a hit at a given position.
       Note that the value returned is a lower bound on the index of that 
       position in the list of hits.
    */
    public int getFirstIndex(int position) {
        if (indexPointers.length == 0) {
            return 1;
        }
        int p = Arrays.binarySearch(indexPositions, position);
        //        System.err.println("First maps " + position + " to " + p);
        if (p >= 0) {
            //            System.err.println("Returning " + indexPointers[p]);
            return indexPointers[p];
        } else {
            p = (p+1)*-1 - 1;
            if (p < 0) { p = 0;}
            //            System.err.println("Returning " + indexPointers[p]);
            return indexPointers[p];
        }
    }
    /* returns the last element in the file of hits that should
       be read when looking for a hit at a given position.
       This is an upper bound on the index of that position in the list of hits.
    */
    public int getLastIndex(int position) {
        if (indexPointers.length == 0) {
            return 1;
        }
        int p = Arrays.binarySearch(indexPositions, position);
        //        System.err.println("Last maps " + position + " to " + p);
        if (p >= 0) {
            p++;
        } else {
            p = (p+1)*-1;
        }
        if (p < indexPointers.length) {
            //            System.err.println("Returning " + indexPointers[p]);
            return indexPointers[p];
        } else {
            //            System.err.println("Returning " + numHits);
            return numHits;
        }
    }
    /* empty constructor for use in saving/restoring from file
     */
    public Header() {}

    /**
     * Creates a header file to index to provided set of hits.
     * The hits must be sorted.
     */
    public Header(IntBuffer hits) {
        init(hits, 16*1024);
    }
    public Header(IntBuffer hits, int pagesize) {
        init(hits, pagesize);
    }
    private void init(IntBuffer hits, int pagesize) {   
        indexPositions = new int[0];
        indexPointers = new int[0];
        numHits = hits.limit();
        int hitstoskip = pagesize / 4;
        ArrayList<Integer> positions = new ArrayList<Integer>();
        ArrayList<Integer> pointers = new ArrayList<Integer>();
        int lastposition = -10;
        int lastpointer = -10 * hitstoskip;
        for (int i = 0; i < numHits; i++) {
            if (i > lastpointer + hitstoskip &&
                lastposition != hits.get(i)) {
                positions.add(hits.get(i));
                pointers.add(i);
                lastpointer = i;
            }
            lastposition = hits.get(i);
        }
        indexPositions = new int[positions.size()];
        indexPointers = new int[positions.size()];
        for (int i = 0; i < indexPositions.length; i++) {
            indexPositions[i] = positions.get(i); 
            indexPointers[i] = pointers.get(i);
        }
    }
    public void writeIndexFile(String fname) throws IOException {
        OutputStream stream = new FileOutputStream(fname);
        byte[] buffer = new byte[8192];
        int[] nh = new int[1];
        nh[0] = numHits;
        Bits.sendInts(nh, stream, buffer);
        Bits.sendInts(indexPositions, stream, buffer);
        Bits.sendInts(indexPointers, stream, buffer);
        stream.close();
    }
    public static Header readIndexFile(String fname) throws IOException {
        File f = new File(fname);
        long size = f.length();
        InputStream stream = new FileInputStream(fname);
        Header h = new Header();
        byte[] buffer = new byte[8192];
        int[] nh = Bits.readInts(1, stream, buffer);
        h.numHits = nh[0];
        size -= 4;
        h.indexPositions = Bits.readInts((int)size/8, stream, buffer);
        h.indexPointers = Bits.readInts((int)size/8, stream, buffer);
        stream.close();

        return h;
    }
 
    public void printIndex() {
        for (int i = 0; i < indexPositions.length; i++) {
            System.err.println(String.format("%d : %d -> %d",i,indexPointers[i],indexPositions[i]));
        }
    }
    
    public void close() {}

}