package edu.mit.csail.cgs.projects.readdb;

import java.nio.*;

/**
 * contains an IntBuffer ib that is
 * derived from the ByteBuffer bb
 */
public class IntBP extends ByteBP {
    IntBuffer ib;
    public IntBP(ByteBuffer b) {
        super(b);
        ib = b.asIntBuffer();
    }
    public IntBP(int size) {
        super(ByteBuffer.allocate(size*4).order(ByteOrder.nativeOrder()));
        ib = bb.asIntBuffer();
    }
    public IntBP slice(int start, int length) {
//         if (length == 0) {
//             System.err.println("LENGTH ZERO");
//         }
        ByteBuffer b;
        synchronized(bb) {
            bb.position(start * 4);
            b = bb.slice();
        }
        b.order(bb.order());
        //        System.err.println("GOING TO LIMIT TO " + (length * 4) + " of " + b.limit());
        b.limit(length * 4);
        IntBP output = new IntBP(b);
        return output;
    }
    public int get(int i) {
        return ib.get(i);
    }
    public void put(int index, int val) {
        ib.put(index,val);
    }
    public int limit() {
        return ib.limit();
    }
    public int size() {
        return ib.limit();
    }
}
