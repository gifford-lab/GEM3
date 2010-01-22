package edu.mit.csail.cgs.projects.readdb;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;

/**
 * Utility methods for reading and writing floats and ints and such
 * in binary form to and from Input and Output Streams.
 *
 */

public class Bits {

    /** Sends the count bytes starting at offset from the ByteBuffer to the channel
     */
    public static void sendBytes(ByteBuffer b, int offset, int count, WritableByteChannel channel) throws IOException {
        synchronized(b) {
            b.position(offset);
            int sent = 0;
            while (sent < count) {
                sent += channel.write(b);
            }
        }
    }
    /** sends the full ByteBuffer to the channel
     */
    public static void sendBytes(ByteBuffer b, WritableByteChannel channel) throws IOException {
        synchronized(b) {
            b.position(0);
            int sent = 0;
            int limit = b.limit();
            int count = 0;
            while (sent < limit) {
                sent += channel.write(b);
                count++;
            }
        }
    }
    /**
     * reads the full ByteBuffer from the channel
     */
    public static void readBytes(ByteBuffer b, ReadableByteChannel channel) throws IOException {
        synchronized (b) {
            b.position(0);
            while (b.remaining() > 0) {
                int r = channel.read(b);
                if (b.remaining() > 0 && r == -1) {
                    throw new IOException("Couldn't fill buffer");
                }
            }
        }
    }
    public static void flipByteOrder(IntBuffer b) {
        for (int i = 0; i < b.limit(); i++) {
            b.put(i, Integer.reverseBytes(b.get(i)));
        }
    }
    public static void flipByteOrder(FloatBuffer b) {
        for (int i = 0; i < b.limit(); i++) {
            b.put(i, Float.intBitsToFloat(Integer.reverseBytes(Float.floatToRawIntBits(b.get(i)))));
        }
    }
    public static ByteBuffer copyByteBuffer(ByteBuffer b) {
        ByteBuffer output = ByteBuffer.allocate(b.limit());
        output.order(b.order());
        for (int i = 0; i < b.limit(); i++) {
            output.put(i,b.get(i));
        }
        return output;
    }


   /** sends the specified integers to the specified output stream.
     *  uses the current value of the order field to determine whether it
     * should send big or little endian.
     *
     * buffer is scratch space
     */
    public static void sendInts(int[] a, OutputStream stream, byte[] buffer, ByteOrder requestedOrder) throws IOException {
        int i = 0;
        ByteBuffer bb = ByteBuffer.wrap(buffer);
        bb.order(requestedOrder);
        while (i < a.length) {
            int end = i + (buffer.length/4) - 1;
            int bufpos = 0;
            if (end >= a.length) {
                end = a.length - 1;
            }            
            for (int j = i; j <= end; j++) {
                bb.putInt(bufpos*4, a[j]);
                bufpos++;
            }		
            stream.write(buffer,0,(end - i + 1) * 4);
            i += buffer.length / 4;
        }
        stream.flush();
    }
    public static void sendFloats(float[] a, OutputStream stream, byte[] buffer, ByteOrder requestedOrder) throws IOException {
        int i = 0;
        ByteBuffer bb = ByteBuffer.wrap(buffer);
        bb.order(requestedOrder);
        while (i < a.length) {
            int end = i + (buffer.length/4) - 1;
            int bufpos = 0;
            if (end >= a.length) {
                end = a.length - 1;
            }            
            for (int j = i; j <= end; j++) {
                bb.putFloat(bufpos*4, a[j]);
                bufpos++;
            }		
            stream.write(buffer,0,(end - i + 1) * 4);
            i += buffer.length / 4;
        }
        stream.flush();
    }
    public static int[] readInts(int count, InputStream instream, byte[] buffer, ByteOrder sendOrder) throws IOException {
        int[] output = new int[count];        
        int outputpos = 0;
        int bytesLeftover = 0;
        ByteBuffer bb = ByteBuffer.wrap(buffer);
        bb.order(sendOrder);
        while (outputpos < count) {
            bb.position(bytesLeftover);
            int toread = Math.min((count - outputpos) * 4, buffer.length - bytesLeftover);
            int bytesavail = instream.read(buffer, bytesLeftover, toread) + bytesLeftover;
            if (bytesavail == -1 && outputpos < count) {
                IOException e = new IOException(String.format("couldn't read enough bytes : %d %d", outputpos, count));
                e.printStackTrace();
                throw e;
            }

            int i = 0;
            for (i = 0; i < bytesavail / 4 && outputpos < count; i++) {
                output[outputpos++] = bb.getInt(i*4);
            }
            int j = i * 4 + 1;
            while (j < bytesavail) {
                buffer[j - i*4] = buffer[j];
                j++;
            }
            bytesLeftover = bytesavail - i*4;
            if (bytesLeftover < 0) {
                System.err.println(String.format("avail was %d but i=%d",bytesavail,i));
            }
            if (bytesLeftover > bb.capacity()) {
                System.err.println(String.format("leftover %d capacity %d", bytesLeftover, bb.capacity()));
            }


        }
        return output;
    }
    public static float[] readFloats(int count, InputStream instream, byte[] buffer, ByteOrder sendOrder) throws IOException {
        float[] output = new float[count];        
        int outputpos = 0;
        int bytesLeftover = 0;
        ByteBuffer bb = ByteBuffer.wrap(buffer);
        bb.order(sendOrder);
        while (outputpos < count) {
            bb.position(bytesLeftover);
            int toread = Math.min((count - outputpos) * 4, buffer.length - bytesLeftover);
            int bytesavail = instream.read(buffer, bytesLeftover, toread) + bytesLeftover;
            int i = 0;
            for (i = 0; i < bytesavail / 4 && outputpos < count; i++) {
                output[outputpos++] = bb.getFloat(i*4);
            }
            int j = i * 4 + 1;
            while (j < bytesavail) {
                buffer[j - i*4] = buffer[j];
                j++;
            }
            bytesLeftover = bytesavail - i*4;
        }
        return output;
    }    
    public static int[] floatToInt(float[] f) {
        int output[] = new int[f.length];
        for (int i = 0; i < f.length; i++) {
            output[i] = Float.floatToRawIntBits(f[i]);
        }
        return output;
    }
    public static float[] intToFloat(int[] i) {
        float[] output = new float[i.length];
        for (int j = 0; j < i.length; j++) {
            output[j] = Float.intBitsToFloat(i[j]);
        }
        return output;
    }


}