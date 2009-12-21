package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.IOException;
import org.junit.*;
import static org.junit.Assert.*;

/** You need to have a server setup with 
 * two accounts.  
 */

public class ReadDBTest {

    private static String hostname, user, passwd, user2, passwd2;
    private static int portnum;

    @Test public void connect() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        c.close();
    }
    @Test(expected=ClientException.class) public void badPassword() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, "foo");
        c.close();
    }
    @Test(expected=ClientException.class) public void badUser() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.close();
    }
    @Test(expected=IOException.class) public void badPort() throws IOException, ClientException {
        Client c = new Client(hostname, portnum+1, user, passwd);
        c.close();
    }
    @Test(expected=ClientException.class) public void badAlignGetChroms() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getChroms("badAlignGetChroms");
        c.close();
    }
    @Test(expected=ClientException.class) public void badAlignGetCount() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getCount("badAlignGetCount","1");
        c.close();
    }
    @Test(expected=ClientException.class) public void badChromGetCount() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getCount("badChromGetCount","1");
        c.close();
    }
    @Test public void storeSorted() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeSorted";
        String chrom = "1";
        int hits[] = {1,2,3};
        float weights[] = {1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hits, weights);
        
        assertTrue(c.exists(name));
        
        Set<String> chroms = c.getChroms(name);
        assertTrue(chroms.size() == 1);
        assertTrue(chroms.contains(chrom));
        assertTrue(String.format("%d != %d", c.getCount(name,chrom), hits.length), c.getCount(name,chrom) == hits.length);
        int read[] = c.getHits(name,chrom);
        assertTrue(read.length == hits.length);
        for (int i = 0; i < read.length; i++) {
            assertTrue(String.format("%d != %d", read[i], hits[i]), read[i] == hits[i]);
        }

        c.close();
    }
    @Test public void storeUnsorted() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeUnsorted";
        String chrom = "1";
        int hits[] = {3,2,1};
        float weights[] = {1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hits,weights);
        
        assertTrue(c.exists(name));
        
        Set<String> chroms = c.getChroms(name);
        assertTrue(chroms.size() == 1);
        assertTrue(chroms.contains(chrom));
        assertTrue(String.format("%d != %d", c.getCount(name,chrom), hits.length), c.getCount(name,chrom) == hits.length);
        int read[] = c.getHits(name,chrom);
        Arrays.sort(hits);
        assertTrue(read.length == hits.length);
        for (int i = 0; i < read.length; i++) {
            assertTrue(read[i] == hits[i]);
        }

        c.close();
    }
    @Test public void storeTwice() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeTwice";
        String chrom = "1";
        int hits[] = {3,2,1};
        float weights[] = {1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hits, weights);
        int hitstwo[] = {6,10,20,30};
        float weightstwo[] = {1.0F, 1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hitstwo, weightstwo);
        int hitsboth[] = {1,2,3,6,10,20,30};
        
        assertTrue(c.exists(name));
        
        Set<String> chroms = c.getChroms(name);
        assertTrue(chroms.size() == 1);
        assertTrue(chroms.contains(chrom));
        if (c.getCount(name,chrom) != hitsboth.length) {
            System.err.println(String.format("Expected %d and got %d", hitsboth.length, c.getCount(name,chrom)));
        }
        assertTrue(String.format("%d != %d", c.getCount(name,chrom), hitsboth.length), c.getCount(name,chrom) == hitsboth.length);
        int read[] = c.getHits(name,chrom);
        Arrays.sort(hitsboth);
        assertTrue(read.length == hitsboth.length);
        for (int i = 0; i < read.length; i++) {
            assertTrue(read[i] == hitsboth[i]);
        }

        c.close();
    }
    @Test public void storeAnotherChrom() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeAnotherChrom";
        String chrom = "1";
        int hits[] = {3,2,1};
        float weights[] = {1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hits, weights);
        int hitstwo[] = {6,10,20};
        c.store(name,"2",hitstwo, weights);
        
        assertTrue(c.exists(name));
        
        Set<String> chroms = c.getChroms(name);
        assertTrue(chroms.size() == 2);
        assertTrue(chroms.contains(chrom));
        assertTrue(chroms.contains("2"));
        assertTrue(String.format("%d != %d", c.getCount(name,chrom), hits.length), c.getCount(name,chrom) == hits.length);
        int read[] = c.getHits(name,chrom);
        Arrays.sort(hits);
        assertTrue(read.length == hits.length);
        for (int i = 0; i < read.length; i++) {
            assertTrue(read[i] == hits[i]);
        }

        c.close();
    }    
    
    @Test public void testACL() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testACL";
        String chrom = "1";
        int hits[] = {3,2,1};
        float weights[] = {1.0F, 1.0F, 1.0F};
        c.store(name,chrom,hits,weights);

        Map<String,Set<String>> acl = c.getACL(name);
        assertTrue(acl.size() == 3);
        assertTrue(acl.get("READ").size() == 1);
        assertTrue(acl.get("WRITE").size() == 1);
        assertTrue(acl.get("ADMIN").size() == 1);
        assertTrue(acl.get("READ").contains(user));
        assertTrue(acl.get("WRITE").contains(user));
        assertTrue(acl.get("ADMIN").contains(user));

        boolean caught = false;
        Client other = null;
        try {
            other = new Client(hostname,portnum,user2,passwd2);
            other.getChroms(name);
        } catch (ClientException e) {
            caught = true;
            other.close();
        }
        assertTrue(caught);

        Set<ACLChangeEntry> changes = new HashSet<ACLChangeEntry>();
        changes.add(new ACLChangeEntry("add","read",user2));
        c.setACL(name, changes);

        acl = c.getACL(name);
        assertTrue(acl.size() == 3);
        assertTrue(acl.get("READ").size() == 2);
        assertTrue(acl.get("WRITE").size() == 1);
        assertTrue(acl.get("ADMIN").size() == 1);
        assertTrue(acl.get("READ").contains(user));
        assertTrue(acl.get("READ").contains(user2));
        assertTrue(acl.get("WRITE").contains(user));
        assertTrue(acl.get("ADMIN").contains(user));

        other = new Client(hostname,portnum,user2,passwd2);
        other.getChroms(name);

        changes.clear();
        changes.add(new ACLChangeEntry("delete","read",user2));
        c.setACL(name, changes);
        acl = c.getACL(name);
        assertTrue(acl.size() == 3);
        assertTrue(acl.get("READ").size() == 1);
        assertTrue(acl.get("WRITE").size() == 1);
        assertTrue(acl.get("ADMIN").size() == 1);
        assertTrue(acl.get("READ").contains(user));
        assertTrue(acl.get("WRITE").contains(user));
        assertTrue(acl.get("ADMIN").contains(user));
        c.close();
    }
    
    @Test public void testRangeQuery() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testRangeQuery";
        String chrom = "1";
        int n = 10000;
        int hits[] = new int[n];
        float weights[] = new float[n];
        for (int i = 0; i < hits.length; i++) {
            hits[i] = (int)Math.round(Math.random() * Integer.MAX_VALUE);
            weights[i] = 1.0F;
        }
        c.store(name,chrom,hits,weights);
        
        Arrays.sort(hits);
        for (int q = 0; q < 30; q++) {
            int start = (int)Math.round(Math.random() * Integer.MAX_VALUE);
            int end = start + (int)(Math.round(Math.random() * Integer.MAX_VALUE) % (Integer.MAX_VALUE - start));
            int count = 0;
            for (int i = 0; i < hits.length; i++) {
                if (hits[i] >= start && hits[i] <= end) {
                    count++;
                }
            }
            System.err.println(String.format("Looking for %d to %d",start,end));
            System.err.println(String.format("My count is %d and server's is %d",count,c.getCountRange(name,chrom,start,end)));
            assertTrue(String.format("%d != %d", count, c.getCountRange(name,chrom,start,end)), count == c.getCountRange(name,chrom,start,end));
            int back[] = c.getHitsRange(name,chrom,start,end);
            int j = 0;
            int k = 0;
            while (j < hits.length && hits[j] < start) { j++;}
            while (j < hits.length && hits[j] <= end) {
                assertTrue(hits[j++] == back[k++]);
            }
        }
        assertTrue(c.getCountRange(name,chrom,0, Integer.MAX_VALUE - 1) == n);
        System.err.println("Done with testRangeQuery");
        c.close();
    }
    @Test public void testHistogram() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testHistogram";
        String chrom = "1";
        int hits[] = {1,1,1,10,10,11,12,25};   
        float weights[] = new float[hits.length];
        for (int i = 0; i < hits.length; i++) { 
            weights[i] = 1.0F;
        }

        System.err.println("Testing histogram");
        c.store(name,chrom,hits,weights);
        System.err.println("Getting histogram");
        
        TreeMap<Integer,Integer> map = c.getHistogram(name,chrom,0,100,10);
        System.err.println(" testing histogram");
        assertTrue(map.size() == 3);
        assertTrue(map.get(5) == 3);
        assertTrue(map.get(15) == 4);
        assertTrue(map.get(25) == 1);
        c.close();
     }
    @Test public void testWeights() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testWeights";
        String chrom = "1";
        int hits[] = {1,1,2,2,3,3,4,4};
        float weights[] = new float[hits.length];
        for (int i = 0; i < hits.length; i++) { 
            if (i % 2 == 0) {
                weights[i] = 0.5F;
            } else {
                weights[i] = 1F;
            }
        }
        c.store(name,chrom,hits,weights);
        int h[] = c.getHits(name,chrom);
        float w[] = c.getWeights(name,chrom);
        assertTrue(w.length == weights.length);
        boolean ok = true;
        for (int i = 0; i < weights.length; i++) {
            if (Math.abs(w[i] - weights[i]) > .00001) {
                System.err.println(String.format("Weight mismatch at %d : %f vs %f (%x vs %x)",
                                                 i, w[i], weights[i], Float.floatToRawIntBits(w[i]), 
                                                 Float.floatToRawIntBits(weights[i])));
                ok = false;
            }
        }
        assertTrue(ok);
        
        h = c.getHitsRange(name,chrom,0,100,0);
        assertTrue(h.length == hits.length);
        h = c.getHitsRange(name,chrom,0,100,.75F);
        assertTrue(h.length == (hits.length / 2));
        h = c.getHitsRange(name,chrom,0,100,.51F);
        assertTrue(h.length == (hits.length / 2));
        h = c.getHitsRange(name,chrom,0,100,.49F);
        assertTrue(h.length == (hits.length));
        c.close();
    }
    @Test public void testHistogram2() throws IOException, ClientException {
        int MAXVALUE = 10000;
        float MAXWEIGHT = 4f;
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testHistogram2";
        for (int q = 0; q < 10; q++) {
            String chrom = Integer.toString(q);
            int[] hits = new int[MAXVALUE];
            float[] weights = new float[MAXVALUE];
            for (int i = 0; i < hits.length; i++) {
                hits[i] = (int)Math.round(Math.random() * MAXVALUE);
                weights[i] = (float)(Math.random() * MAXWEIGHT);
            }
            Arrays.sort(hits);
            c.store(name,chrom,hits,weights);           

            for (int q2 = 0; q2 < 10; q2++) {
                int start = (int)Math.round(Math.random() * MAXVALUE);
                int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
                float weight = (float)Math.random() * MAXWEIGHT;
                int binsize = 10 + (int)Math.round(Math.random() * 40);
                if (end < start) {
                    throw new RuntimeException("end < start");
                }
                Map<Integer,Integer> plainhist = c.getHistogram(name,chrom,start,end,binsize);
                Map<Integer,Integer> minweighthist = c.getHistogram(name,chrom,start,end,binsize,weight,0);
                Map<Integer,Float> weighthist = c.getWeightHistogram(name,chrom,start,end,binsize);
                Map<Integer,Float> minwwhist = c.getWeightHistogram(name,chrom,start,end,binsize,weight,0);

                Map<Integer,Integer> myplain = new HashMap<Integer,Integer>();
                Map<Integer,Integer> myminw = new HashMap<Integer,Integer>();
                Map<Integer,Float> myw = new HashMap<Integer,Float>();
                Map<Integer,Float> myww = new HashMap<Integer,Float>();
                for (int i = 0; i < hits.length; i++) {
                    int bin = ((hits[i] - start) / binsize) * binsize + binsize / 2 + start;
                    if (hits[i] >= start && hits[i] <= end) {
                        if (myplain.containsKey(bin)) {
                            myplain.put(bin,1 + myplain.get(bin));
                        } else {
                            myplain.put(bin,1);
                        }
                        if (myw.containsKey(bin)) {
                            myw.put(bin,weights[i] + myw.get(bin));
                        } else {
                            myw.put(bin,weights[i]);
                        }                        
                        if (weights[i] >= weight) {
                            if (myminw.containsKey(bin)) {
                                myminw.put(bin,1 + myminw.get(bin));
                            } else {
                                myminw.put(bin,1);
                            }
                            if (myww.containsKey(bin)) {
                                myww.put(bin,weights[i] + myww.get(bin));
                            } else {
                                myww.put(bin,weights[i]);
                            }                            
                        }
                    }
                }
//                 System.err.println(String.format("Start %d, Stop %d, Binsize %d", start,end,binsize));
//                 System.err.println("My plain keys " + myplain.keySet());
//                 System.err.println("Plain keys " + plainhist.keySet());
                assertTrue(myplain.size() == plainhist.size());
                assertTrue(myminw.size() == minweighthist.size());
                assertTrue(myw.size() == weighthist.size());
                assertTrue(myww.size() == minwwhist.size());
                for (int bin : myplain.keySet()) {
                    assertTrue(String.format("plain bin %d : %d != %d", bin, myplain.get(bin), plainhist.get(bin)), 
                               myplain.get(bin) == plainhist.get(bin));
                }
                for (int bin : myminw.keySet()) {
                    assertTrue(String.format("minw bin %d : %d != %d", bin, myminw.get(bin), minweighthist.get(bin)), 
                               myminw.get(bin) == minweighthist.get(bin));                   
                }
                for (int bin : myw.keySet()) {
                    assertTrue(String.format("weight bin %d : %f != %f", bin, myw.get(bin), weighthist.get(bin)), 
                               Math.abs(myw.get(bin) - weighthist.get(bin)) < .0001);
                }
                for (int bin : myww.keySet()) {
                    assertTrue(String.format("ww bin %d : %f != %f", bin, myww.get(bin), minwwhist.get(bin)), 
                               Math.abs(myww.get(bin) - minwwhist.get(bin)) < .0001);
                }
            }
        }        
        c.close();
    }
    


    public static void main(String args[]) {
        hostname = args[0];
        portnum = Integer.parseInt(args[1]);
        user = args[2];
        passwd = args[3];
        user2 = args[4];
        passwd2 = args[5];

        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.ReadDBTest");
    }
}