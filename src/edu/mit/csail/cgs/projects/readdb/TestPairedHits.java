package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.*;
import org.junit.*;
import static org.junit.Assert.*;

public class TestPairedHits {

    private static int chrom;
    private static String prefix;
    private static boolean isLeft;
    private Header header;
    private PairedHits hitsfile;
    private List<PairedHit> hits;

    public TestPairedHits() throws IOException {
        int len = 1000;
        chrom = (int)(Math.random() * 10000);
        hits = new ArrayList<PairedHit>();
        for (int i = 0; i < len; i++) {
            hits.add(new PairedHit(isLeft ? chrom : (int)(Math.random() * 1000),
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   isLeft ? (int)(Math.random() * 1000) : chrom,
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   (float)Math.random()));
        }        
        Collections.sort(hits, isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator());
        PairedHits.writePairedHits(hits, prefix, chrom, isLeft);
        
        hitsfile = new PairedHits(prefix,chrom,isLeft);
        header = new Header(hitsfile.getPositionsBuffer().ib);
        header.writeIndexFile(prefix + chrom + ".index");        
        header = Header.readIndexFile(prefix + chrom + ".index");
    }

    @Test public void testCount() {
        assertEquals(hits.size(), header.getNumHits());
    }

    @Test public void testArraySort() {
        List<PairedHit> ph = new ArrayList<PairedHit>();
        for (int i = 0; i < 1000; i++) {
            ph.add(new PairedHit(chrom,
                                 (int)(Math.random() * 100000000),
                                 i % 3 == 1,
                                 (short)(Math.random()*100+1),
                                 chrom,
                                 (int)(Math.random() * 100000000),
                                 i % 3 == 1,
                                 (short)(Math.random()*100+1),
                                 (float)Math.random()));
        }        
        Collections.sort(ph, new PairedHitLeftComparator());
        for (int i = 1; i < hits.size(); i++) {
            assertTrue(ph.get(i-1).leftPos <= ph.get(i).leftPos);
        }
        Collections.sort(ph, new PairedHitRightComparator());
        for (int i = 1; i < hits.size(); i++) {
            assertTrue(ph.get(i-1).rightPos <= ph.get(i).rightPos);
        }        
    }

    @Test public void testFileSort() {
        IntBP positions = hitsfile.getPositionsBuffer();
        for (int i = 1; i < hits.size(); i++) {
            assertTrue(positions.get(i-1) < positions.get(i));
        }        
    }

    @Test public void testChroms() {
        IntBP chroms = hitsfile.getChromsBuffer();
        for (int i = 0; i < hits.size(); i++) {
            assertEquals(hitsfile.isLeft() ? hits.get(i).rightChrom : hits.get(i).leftChrom, chroms.get(i));
        }
    }

    @Test public void testPositions() {
        IntBP positions = hitsfile.getOtherPosBuffer();
        for (int i = 0; i < hits.size(); i++) {
            assertEquals(hitsfile.isLeft() ? hits.get(i).rightPos : hits.get(i).leftPos, positions.get(i));
        }
    }

    @Test public void testStrand() {
        IntBP las = hitsfile.getLASBuffer();
        for (int i = 0; i < hits.size(); i++) {
            assertEquals(!hitsfile.isLeft() ? hits.get(i).rightStrand : hits.get(i).leftStrand, Hits.getStrandOne(las.get(i)));
            assertEquals(hitsfile.isLeft() ? hits.get(i).rightStrand : hits.get(i).leftStrand, Hits.getStrandTwo(las.get(i)));
        }
    }
    @Test public void testLength() {
        IntBP las = hitsfile.getLASBuffer();
        for (int i = 0; i < hits.size(); i++) {
            assertEquals(!hitsfile.isLeft() ? hits.get(i).rightLength : hits.get(i).leftLength, Hits.getLengthOne(las.get(i)));
            assertEquals(hitsfile.isLeft() ? hits.get(i).rightLength : hits.get(i).leftLength, Hits.getLengthTwo(las.get(i)));
        }
    }
    @Test public void testAppend () throws IOException {
        List<PairedHit> orighits = new ArrayList<PairedHit>();
        List<PairedHit> newhits = new ArrayList<PairedHit>();
        int len = 1000;
        for (int i = 0; i < len; i++) {
            orighits.add(new PairedHit(isLeft ? chrom : (int)(Math.random() * 1000),
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   isLeft ? (int)(Math.random() * 1000) : chrom,
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   (float)Math.random()));
        }        
        for (int i = 0; i < len; i++) {
            newhits.add(new PairedHit(isLeft ? chrom : (int)(Math.random() * 1000),
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   isLeft ? (int)(Math.random() * 1000) : chrom,
                                   (int)(Math.random() * 100000000),
                                   i % 3 == 1,
                                   (short)(Math.random()*100+1),
                                   (float)Math.random()));
        }      
        Collections.sort(orighits, isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator());
        Collections.sort(newhits, isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator());
        PairedHits.writePairedHits(orighits, prefix, chrom, isLeft);
        hitsfile = new PairedHits(prefix,chrom,isLeft);
        hitsfile.appendPairedHits(newhits,prefix + "appended", chrom, isLeft);
        hitsfile = new PairedHits(prefix + "appended",chrom,isLeft);

        assertEquals(hitsfile.isLeft(), isLeft);
        ArrayList<PairedHit> allhits = new ArrayList<PairedHit>();
        allhits.addAll(orighits);
        allhits.addAll(newhits);
        Collections.sort(allhits, isLeft ? new PairedHitLeftComparator() : new PairedHitRightComparator());

        IntBP positions = hitsfile.getPositionsBuffer();
        IntBP las = hitsfile.getLASBuffer();
        FloatBP weights = hitsfile.getWeightsBuffer();
        IntBP otherchrom = hitsfile.getChromsBuffer();
        IntBP otherpos = hitsfile.getOtherPosBuffer();

        for (int i = 0; i < allhits.size(); i++) {
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)), 
                         allhits.get(i).leftChrom, isLeft ? chrom : otherchrom.get(i));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).leftPos, isLeft ? positions.get(i) : otherpos.get(i));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).leftStrand, isLeft ? Hits.getStrandOne(las.get(i)) : Hits.getStrandTwo(las.get(i)));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).leftLength, isLeft ? Hits.getLengthOne(las.get(i)) : Hits.getLengthTwo(las.get(i)));

            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).rightChrom, !isLeft ? chrom : otherchrom.get(i));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).rightPos, !isLeft ? positions.get(i) : otherpos.get(i));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).rightStrand, !isLeft ? Hits.getStrandOne(las.get(i)) : Hits.getStrandTwo(las.get(i)));
            assertEquals(String.format("Hit %d, %s",i,allhits.get(i)),
                         allhits.get(i).rightLength, !isLeft ? Hits.getLengthOne(las.get(i)) : Hits.getLengthTwo(las.get(i)));
            
            assertTrue(String.format("Hit %d, %s",i,allhits.get(i)),
                       Math.abs(allhits.get(i).weight - weights.get(i)) < .00001);                               
        }


    }

   public static void main(String args[]) {
        prefix = args[0];
        if (!prefix.endsWith(System.getProperty("file.separator"))) {
            prefix = prefix + System.getProperty("file.separator");
        }
        chrom = Integer.parseInt(args[1]);
        isLeft = Boolean.parseBoolean(args[2]);

        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.TestPairedHits");
    }
}