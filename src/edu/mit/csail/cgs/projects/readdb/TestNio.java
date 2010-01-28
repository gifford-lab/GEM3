package edu.mit.csail.cgs.projects.readdb;

import java.nio.*;


public class TestNio {


    public static void main(String args[]) {
        ByteBuffer bb = ByteBuffer.allocate(4);
        bb.order(ByteOrder.nativeOrder());
        System.out.println("Native order is " + bb.order());
        IntBuffer ib = bb.asIntBuffer();
        ib.put(1);
        System.out.println(String.format("%d %d %d %d", bb.get(0), bb.get(1), bb.get(2), bb.get(3)));
        System.err.println("Byte was " + ib.get(0));
        System.err.println(ib.order() + "   " + bb.order());

    }    

}