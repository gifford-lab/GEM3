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
        int[] temp = new int[NUMHITS];
        for (int i = 0; i < NUMHITS; i++) {
            temp[i] =  (int)Math.round(Math.random() * 4000);
        }
        Arrays.sort(temp);
        for (int i = 0; i < NUMHITS; i++) {
            hits.put(i, temp[i]);
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

    @Test public void testLAS() {
        for (short s = 0; s < 1000; s += 10) {
            int l = Hits.makeLAS(s, true);
            assertTrue(Hits.getStrandOne(l));
            assertEquals(Hits.getLengthOne(l), s);

            l = Hits.makeLAS(s, false);
            assertFalse(Hits.getStrandOne(l));
            assertEquals(Hits.getLengthOne(l), s);
        }
        for (short s = 0; s < 1000; s += 10) {
            for (short t = 0; t < 1000; t += 10) {
                int l = Hits.makeLAS(s,true,t,true);
                assertTrue(Hits.getStrandOne(l));
                assertEquals(Hits.getLengthOne(l), s);
                assertTrue(Hits.getStrandTwo(l));
                assertEquals(Hits.getLengthTwo(l), t);

                l = Hits.makeLAS(s,true,t,false);
                assertTrue(Hits.getStrandOne(l));
                assertEquals(Hits.getLengthOne(l), s);
                assertFalse(Hits.getStrandTwo(l));
                assertEquals(Hits.getLengthTwo(l), t);

                l = Hits.makeLAS(s,false,t,true);
                assertFalse(Hits.getStrandOne(l));
                assertEquals(Hits.getLengthOne(l), s);
                assertTrue(Hits.getStrandTwo(l));
                assertEquals(Hits.getLengthTwo(l), t);

                l = Hits.makeLAS(s,false,t,false);
                assertFalse(Hits.getStrandOne(l));
                assertEquals(Hits.getLengthOne(l), s);
                assertFalse(Hits.getStrandTwo(l));
                assertEquals(Hits.getLengthTwo(l), t);
            }
        }               
    }

    @Test public void testGetAllHits (){
        IntBP positions = hitsfile.getPositionsBuffer();
        FloatBP w = hitsfile.getWeightsBuffer();
        IntBP l = hitsfile.getLASBuffer();
        boolean fail = false;
        for (int i = 0; i < hits.size(); i++) {            
            if (positions.get(i) != hits.get(i)) {
                System.err.println(String.format("At %d, %d != %d", i, positions.get(i), hits.get(i)));
                fail = true;
            }
            if (Math.abs(w.get(i) - weights.get(i)) > .001) {
                System.err.println(String.format("At %d, %.4f != %.4f",i,weights.get(i), w.get(i)));
                fail = true;
            }
            if (l.get(i) != las.get(i)) {
                System.err.println(String.format("At %d, las %d != %d", i, las.get(i), l.get(i)));
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
        assertTrue(String.format("cornercase : %d,%d, %d,%d -> %d,%d",
                                 start,end, firstindex, lastindex, p[0], p[1]), p[1] == header.getNumHits());

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
                                                 start,end,binsize,0,null,null,false);
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

    @Test public void testDeDupHistogram() throws IOException {
        for (int q = 0; q < 300; q++) {
            int start = (int)Math.round(Math.random() * (MAXVALUE - 10));
            int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
            if (end < start) {
                throw new RuntimeException("end < start");
            }
            int binsize = 10 + (int)Math.round(Math.random() * 40);
            int k = 0;
            // while (hits.get(k) < binsize) {
            //     System.out.println("k=" + k + "  val=" + hits.get(k));
            //     k++;
            // }


            int[] histogram = hitsfile.histogram(header.getFirstIndex(start),
                                                 header.getLastIndex(end),
                                                 start,end,binsize,1,null,null,false);

            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end) {
                    if (i == 0 || hits.get(i-1) != hits.get(i)) {
                        myhist[(hits.get(i) - start) / binsize]++;
                    }
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                if (myhist[i] != histogram[i]) {
                    System.err.println(String.format("At dup1  %d, %d != %d",i,myhist[i],histogram[i]));
                }
                assertTrue(myhist[i] == histogram[i]);
            }


            histogram = hitsfile.histogram(header.getFirstIndex(start),
                                                 header.getLastIndex(end),
                                                 start,end,binsize,2,null,null,false);

            myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end) {
                    if (i < 2  || (hits.get(i-1) != hits.get(i) ||
                                   hits.get(i-2) != hits.get(i))) {
                        myhist[(hits.get(i) - start) / binsize]++;
                    }
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                if (myhist[i] != histogram[i]) {
                    System.err.println(String.format("At dup2  %d, %d != %d",i,myhist[i],histogram[i]));
                }
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
                                                 start,end,binsize,0,weight,null,false);
            int[] myhist = new int[(end - start) / binsize + 1];
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i) >= start && hits.get(i) <= end && weights.get(i) >= weight) {
                    myhist[(hits.get(i) - start) / binsize]++;
                }
            }
            assertTrue(myhist.length == histogram.length);
            for (int i = 0; i < myhist.length; i++) {
                assertTrue(Math.abs(myhist[i] - histogram[i]) < .001);
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
                                                         start,end,binsize,0,null,null,false);
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
        if (!prefix.endsWith(System.getProperty("file.separator"))) {
            prefix = prefix + System.getProperty("file.separator");
        }

        chrom = Integer.parseInt(args[1]);
        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.TestHits");
    }

}