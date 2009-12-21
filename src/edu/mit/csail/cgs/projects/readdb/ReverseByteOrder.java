package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * Reverse the byte-order in a file of floats or ints.  You only
 * need this if you're Alex and the old version of your code 
 * wrote the files with the wrong byte order
 */

public class ReverseByteOrder {

    public static void main(String fnames[]) {
        for (int i = 0; i < fnames.length; i++) {
            String fname = fnames[i];
            try {
                RandomAccessFile raf = new RandomAccessFile(fname, "rw");
                FileChannel fc = raf.getChannel();
                ByteBuffer bb = fc.map(FileChannel.MapMode.READ_WRITE,
                                       0,
                                       fc.size());            
                IntBuffer ib = bb.asIntBuffer();
                Bits.flipByteOrder(ib);
                ib = null;
                bb = null;
                fc.close();
                raf.close();
                System.out.println("fixed " + fname);
            } catch (Exception e) {
                System.out.println("failed " + fname);
                e.printStackTrace();
            }
        }
    }

}