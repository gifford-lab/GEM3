package edu.mit.csail.cgs.projects.readdb;

import java.nio.*;

/**
 * Ancestor class for wrappers that combine
 * a ByteBuffer and a view of the ByteBuffer such
 * as an IntBuffer.
 *
 * We need this combination because some methods are most convenient on
 * the view (eg, IntBuffer.get() is much better than ByteBuffer.getInt
 * because it reduces the chance of making mistakes with the index)
 * but sometimes you need the ByteBuffer too, for example, when sending
 * across the network.
 *
 * BP = buffer pair
 */
public class ByteBP {

    protected ByteBuffer bb;
    public ByteBP(ByteBuffer b) {
        bb = b;
    }

}