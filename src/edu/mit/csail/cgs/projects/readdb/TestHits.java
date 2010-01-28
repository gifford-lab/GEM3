package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.*;
import java.nio.*;
import org.junit.*;
import static org.junit.Assert.*;


public class TestHits {

    private IntBP hits, las;
    private FloatBP weights;
    private Header header;
    private SingleHits hitsfile;
    private static String prefix;
    private static int chrom;
    private static int MAXVALUE;
    private static float MAXWEIGHT = 4f;
    private static int NUMHITS = 100000;

    public TestHits() throws IOException {
        hits = new IntBP(NUMHITS);
        weights = new FloatBP(NUMHITS);
        las = new IntBP(NUMHITS);
        hits.put(0, (int)Math.round(Math.random() * 5));
        for (int i = 1; i < NUMHITS; i++) {
            hits.put(i, hits.get(i-1) + (int)Math.round(Math.random() * 4));
        }
        MAXVALUE = hits.get(NUMHITS - 1);
        for (int i = 0; i < hits.limit(); i++) {
            weights.put(i,(float)(Math.random() * MAXWEIGHT));
            las.put(i, Hits.makeLAS((short)(hits.get(i) % 500), hits.get(i) % 2 == 1,
                                    (short)(hits.get(i) % 1000), hits.get(i) % 5 == 1));
        }
        SingleHits.writeSingleHits(hits, weights, las, prefix, chrom);
        hitsfile = new SingleHits(prefix, chrom);
        header = new Header(hitsfile.getPositionsBuffer().ib);
        header.writeIndexFile(prefix + chrom + ".index");

        header = Header.readIndexFile(prefix + chrom + ".index");
    }

    @Test public void testCount() {
        assertTrue(hits.size() == header.getNumHits());        
    }

    @Test public void testGetAllHits (){
        IntBP positions = hitsfile.getPositionsBuffer();
        boolean fail = false;
        for (int i = 0; i < hits.size(); i++) {            
            if (positions.get(i) != hits.get(i)) {
                System.err.println(String.format("At %d, %d != %d", i, positions.get(i), hits.get(i)));
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

        start = hits.get(1000);
        firstindex = header.getFirstIndex(start);
        lastindex = header.getLastIndex(end);
        p = hitsfile.getIndices(firstindex,lastindex,start,end);
        assertTrue(p[0] <= p[1]);        
        assertTrue(p[1] == header.getNumHits());

        start = -10;
        end = hits.get(1000);
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
            assertTrue(p[1] <= hits.size());
            boolean ok = true;
            if (p[0] == p[1]) {
                for (int i = 0; i < hits.size(); i++) {
                    if (hits.get(i) >= start && hits.get(i) <= end) {
                        ok = false;
                    }
                }
            } else {
                ok = (hits.get(p[0]) >= start && (p[0] == 0 || hits.get(p[0]-1) < start)) && 
                    (p[1] == hits.size() || hits.get(p[1]) > end && hits.get(p[1]-1) <= end);
            }
            if (!ok) {
                System.err.println(String.format("Not OK : getIndices(%d,%d,%d,%d) -> %d,%d",
                                                 firstindex,lastindex,start,end,p[0],p[1]));
                System.err.println(String.format("      %d,%d,%d .. %d, %d, %d",
                                                 (p[0] == 0 ? -10 : hits.get(p[0]-1)),hits.get(p[0]), hits.get(p[0]+1), hits.get(p[1]-1), (p[1] == hits.size() ? -10 : hits.get(p[1])), (p[1] == hits.size() - 1 ? -10 : hits.get(p[1]+1))));
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
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end) {
                    count++;
                }
            }
            int reported = hitsfile.getCountBetween(header.getFirstIndex(start),
                                                    header.getLastIndex(end),
                                                    start,end, null, null);
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
                                                     start,end, null, null);
            FloatBP wreported = hitsfile.getWeightsBetween(header.getFirstIndex(start),
                                                           header.getLastIndex(end),
                                                           start,end, null, null);
            int k = 0; int j = 0;
            while (k < hits.size() && hits.get(k) < start) { k++;}
            while (k < hits.size() && hits.get(k) <= end) {
                if (hits.get(k) != reported.get(j)) {
                    System.err.println(String.format("Mismatch %d (%d) != %d (%d)",
                                                     hits.get(k), k, reported.get(j), j));
                }
                assertTrue(hits.get(k) == reported.get(j));
                assertTrue(weights.get(k) == wreported.get(j));
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
                                                 start,end,binsize,null,null,false);
            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end) {
                    myhist[(hits.get(i) - start) / binsize]++;
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
                                                 start,end,binsize,weight,null,false);
            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end && weights.get(i) >= weight) {
                    myhist[(hits.get(i) - start) / binsize]++;
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
                                                         start,end,binsize,null,null,false);
            float[] myhist = new float[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end) {
                    myhist[(hits.get(i) - start) / binsize] += weights.get(i);
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                assertTrue(myhist[i] == histogram[i]);
            }
        }                
    }


    public static void main(String args[]) {
        prefix = args[0];
        chrom = Integer.parseInt(args[1]);
        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.TestHits");
    }

}