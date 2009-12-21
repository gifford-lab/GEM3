package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.*;
import java.nio.*;
import org.junit.*;
import static org.junit.Assert.*;


public class HitsTest {

    private int hits[];
    private float weights[];
    private Header header;
    private Hits hitsfile;
    private static String hitsFilename, weightsFilename;
    private static int MAXVALUE = 10000;
    private static float MAXWEIGHT = 4f;
    private static int NUMHITS = 100000;

    public HitsTest() throws IOException {
        hits = new int[NUMHITS];
        weights = new float[NUMHITS];
        for (int i = 0; i < hits.length; i++) {
            hits[i] = (int)Math.round(Math.random() * MAXVALUE);
            weights[i] = (float)(Math.random() * MAXWEIGHT);
        }
        Arrays.sort(hits);
        Bits.sendInts(hits, 
                      new FileOutputStream(hitsFilename),
                      new byte[8192],
                      ByteOrder.nativeOrder());
        Bits.sendFloats(weights, 
                        new FileOutputStream(weightsFilename),
                        new byte[8192],
                        ByteOrder.nativeOrder());
        header = new Header(IntBuffer.wrap(hits), 40);
        hitsfile = new Hits(hitsFilename, weightsFilename);
    }

    @Test public void testCount() {
        assertTrue(hits.length == header.getNumHits());        
    }

    @Test public void testGetAllHits (){
        IntBP positions = hitsfile.getPositionsBuffer();
        boolean fail = false;
        for (int i = 0; i < hits.length; i++) {            
            if (positions.get(i) != hits[i]) {
                System.err.println(String.format("At %d, %d != %d", i, positions.get(i), hits[i]));
                fail = true;
            }
        }
        assertTrue(!fail);
    }

    @Test public void testCornerCases() {
        int start = -100;
        int end = -20;
        int firstindex = header.getFirstIndex(start);
        int lastindex = header.getLastIndex(end);
        int p[] = hitsfile.getIndices(firstindex,lastindex,start,end);
        assertTrue(p[0] >= p[1]);        

        start = MAXVALUE+10;
        end = MAXVALUE+10;
        firstindex = header.getFirstIndex(start);
        lastindex = header.getLastIndex(end);
        p = hitsfile.getIndices(firstindex,lastindex,start,end);
        assertTrue(p[0] >= p[1]);        

        start = hits[1000];
        firstindex = header.getFirstIndex(start);
        lastindex = header.getLastIndex(end);
        p = hitsfile.getIndices(firstindex,lastindex,start,end);
        assertTrue(p[0] <= p[1]);        
        assertTrue(p[1] == header.getNumHits());

        start = -10;
        end = hits[1000];
        firstindex = header.getFirstIndex(start);
        lastindex = header.getLastIndex(end);
        p = hitsfile.getIndices(firstindex,lastindex,start,end);
        assertTrue(p[0] == 0);
    }

    @Test public void testGetIndices() {
        boolean allok = true;
        for (int q = 1; q < 5000; q++) {
            int start = (int)Math.round(Math.random() * MAXVALUE);
            int end = start + (int)Math.round(Math.random() * (MAXVALUE - start));
            int firstindex = header.getFirstIndex(start);
            int lastindex = header.getLastIndex(end);
            
            int p[] = hitsfile.getIndices(firstindex,lastindex,start,end);
            assertTrue(p[0] >= 0);
            assertTrue(p[1] <= hits.length);
            boolean ok = true;
            if (p[0] == p[1]) {
                for (int i = 0; i < hits.length; i++) {
                    if (hits[i] >= start && hits[i] <= end) {
                        ok = false;
                    }
                }
            } else {
                ok = (hits[p[0]] >= start && (p[0] == 0 || hits[p[0]-1] < start)) && 
                    (p[1] == hits.length || hits[p[1]] > end && hits[p[1]-1] <= end);
            }
            if (!ok) {
                System.err.println(String.format("Not OK : getIndices(%d,%d,%d,%d) -> %d,%d",
                                                 firstindex,lastindex,start,end,p[0],p[1]));
                System.err.println(String.format("      %d,%d,%d .. %d, %d, %d",
                                                 (p[0] == 0 ? -10 : hits[p[0]-1]),hits[p[0]], hits[p[0]+1], hits[p[1]-1], (p[1] == hits.length ? -10 : hits[p[1]]), (p[1] == hits.length - 1 ? -10 : hits[p[1]+1])));
            }
            allok = allok && ok;
        }
        assertTrue(allok);
    }
    
    @Test public void testCountRange() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * MAXVALUE);
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }

            int count = 0;
            for (int i = 0; i < hits.length; i++) {
                if (hits[i] >= start && hits[i] <= end) {
                    count++;
                }
            }
            int reported = hitsfile.getCountBetween(header.getFirstIndex(start),
                                                    header.getLastIndex(end),
                                                    start,end);
            assertTrue(count == reported);
        }
    }
    
    @Test public void testGetHits() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * MAXVALUE);
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }
            IntBP reported = hitsfile.getHitsBetween(header.getFirstIndex(start),
                                                     header.getLastIndex(end),
                                                     start,end);
            FloatBP wreported = hitsfile.getWeightsBetween(header.getFirstIndex(start),
                                                           header.getLastIndex(end),
                                                           start,end);
            int k = 0; int j = 0;
            while (k < hits.length && hits[k] < start) { k++;}
            while (k < hits.length && hits[k] <= end) {
                if (hits[k] != reported.get(j)) {
                    System.err.println(String.format("Mismatch %d (%d) != %d (%d)",
                                                     hits[k], k, reported.get(j), j));
                }
                assertTrue(hits[k] == reported.get(j));
                assertTrue(weights[k] == wreported.get(j));
                k++;
                j++;
            }
        }
    }

    @Test public void testHistogram() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * (MAXVALUE - 10));
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }
            int binsize = 10 + (int)Math.round(Math.random() * 40);
            int[] histogram = hitsfile.histogram(header.getFirstIndex(start),
                                                 header.getLastIndex(end),
                                                 start,end,binsize,0);
            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.length; i++) {
                if (hits[i] >= start && hits[i] <= end) {
                    myhist[(hits[i] - start) / binsize]++;
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                assertTrue(myhist[i] == histogram[i]);
            }
        }        
    }

    @Test public void testMinWeightHistogram() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * MAXVALUE);
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }
            int binsize = 10 + (int)Math.round(Math.random() * 40);
            float weight = (float)Math.random() * MAXWEIGHT;
            int[] histogram = hitsfile.histogram(header.getFirstIndex(start),
                                                 header.getLastIndex(end),
                                                 start,end,binsize,weight,0);
            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.length; i++) {
                if (hits[i] >= start && hits[i] <= end && weights[i] >= weight) {
                    myhist[(hits[i] - start) / binsize]++;
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                assertTrue(myhist[i] == histogram[i]);
            }
        }        
    }
    
    @Test public void testWeightHistogram() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * MAXVALUE);
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }
            int binsize = 10 + (int)Math.round(Math.random() * 40);
            float[] histogram = hitsfile.weightHistogram(header.getFirstIndex(start),
                                                         header.getLastIndex(end),
                                                         start,end,binsize,0);
            float[] myhist = new float[(end - start) / binsize + 1];
            for (int i = 0; i < hits.length; i++) {
                if (hits[i] >= start && hits[i] <= end) {
                    myhist[(hits[i] - start) / binsize] += weights[i];
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                assertTrue(myhist[i] == histogram[i]);
            }
        }                
    }


    public static void main(String args[]) {
        hitsFilename = args[0];
        weightsFilename = args[1];
        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.HitsTest");
    }

}