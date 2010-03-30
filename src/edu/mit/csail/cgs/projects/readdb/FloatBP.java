package edu.mit.csail.cgs.projects.readdb;

import java.nio.*;

/**
 * contains a FloatBuffer ib that is
 * derived from the ByteBuffer bb
 */
public class FloatBP extends ByteBP {
    protected FloatBuffer fb;
    public FloatBP(ByteBuffer b) {
        super(b);
        fb = b.asFloatBuffer();
    }
    public FloatBP(int size) {
        super(ByteBuffer.allocate(size*4));
        fb = bb.asFloatBuffer();
    }
    public FloatBP slice(int start, int length) {
        ByteBuffer b;
        synchronized(bb) {
            bb.position(start * 4);
            b = bb.slice();
        }
        b.limit(length * 4);
        return new FloatBP(b);
    }
    public float get(int i) {
        return fb.get(i);
    }
    public void put(int index, float val) {
        fb.put(index,val);
    }
    public int limit() {
        return fb.limit();
    }
    public int size() {
        return fb.limit();
    }
}
