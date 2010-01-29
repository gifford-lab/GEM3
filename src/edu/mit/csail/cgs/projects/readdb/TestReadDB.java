package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.IOException;
import org.junit.*;
import static org.junit.Assert.*;

/**
 * Regression tests for readb.  Uses the Client class to communicate with a server
 * (should be newly initialized)
 * You need to have a server setup with two accounts. 
 */

public class TestReadDB {

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
    @Test(expected=ClientException.class) public void badAlignGetChromsPaired() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getChroms("badAlignGetChroms",true,true);
        c.close();
    }
    @Test(expected=ClientException.class) public void badAlignGetChromsSingle() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getChroms("badAlignGetChroms",false,true);
        c.close();
    }
    @Test(expected=ClientException.class) public void badAlignGetCountSingle() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getCount("badAlignGetCount",false,true,true);
        c.close();
    }
    @Test(expected=ClientException.class) public void badAlignGetCountPaired() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, "foo", passwd);
        c.getCount("badAlignGetCount",true,true,true);
        c.close();
    }
    @Test public void storeSortedSingle() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeSortedSingle";
        ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
        hits.add(new SingleHit(1, 1, 1.0F, true, 10));
        hits.add(new SingleHit(1, 2, 2.0F, true, 20));
        hits.add(new SingleHit(1, 3, 3.0F, true, 30));
        hits.add(new SingleHit(1, 1000, 4.0F, true, 10));
        hits.add(new SingleHit(10, 10, 4.0F, true, 10));
        hits.add(new SingleHit(10, 20, 4.0F, true, 10));
        hits.add(new SingleHit(10, 40, 4.0F, false, 10));

        c.storeSingle(name,hits);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertTrue(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(), 2);
        
        assertEquals(name + " chrom 1 has 4 hits (either strand)", c.getCount(name, 1, false,null,null,null,null,null), 4);
        assertEquals(name + " chrom 10 has 4 hits (either strand)", c.getCount(name, 10, false,null,null,null,null,null), 3);
        assertEquals(name + " chrom 1 has 4 hits (+ strand)", c.getCount(name, 1, false,null,null,null,null,true), 4);
        assertEquals(name + " chrom 10 has 4 hits (+ strand)", c.getCount(name, 10, false,null,null,null,null,true), 2);
        assertEquals(name + " chrom 1 has 0 hits (- strand)", c.getCount(name, 1, false,null,null,null,null,false), 0);
        assertEquals(name + " chrom 10 has 1 hits (- strand)", c.getCount(name, 10, false,null,null,null,null,false), 1);

        List<SingleHit> retrieved = c.getSingleHits(name,1,null,null,null,null);
        assertEquals(name + " get chr 1 hits", retrieved.size(), 4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,true);
        assertEquals(name + " get chr 1 hits +", retrieved.size(), 4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits + pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,false);
        assertEquals(name + " get chr 1 hits -", retrieved.size(), 0);

        retrieved = c.getSingleHits(name,10,null,null,null,null);
        assertEquals(name + " get chr 10 hits", retrieved.size(), 3);
        for (int i = 0; i < 3; i++) {
            assertTrue(name + " get chr 10 hits pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,true);
        assertEquals(name + " get chr 10 hits +", retrieved.size(), 2);
        for (int i = 0; i < 2; i++) {
            assertTrue(name + " get chr 10 hits + pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,false);
        assertEquals(name + " get chr 10 hits -", retrieved.size(), 1);
        for (int i = 0; i < 1; i++) {
            assertTrue(name + " get chr 10 hits - pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+6)));
        }

        c.close();
    }
    @Test public void storeTwiceSingle() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeTwiceSingle";
        ArrayList<SingleHit> hitsOne = new ArrayList<SingleHit>();
        ArrayList<SingleHit> hitsTen = new ArrayList<SingleHit>();
        hitsOne.add(new SingleHit(1, 1, 1.0F, true, 10));
        hitsOne.add(new SingleHit(1, 1000, 4.0F, true, 10));
        hitsOne.add(new SingleHit(1, 3, 3.0F, true, 30));
        hitsOne.add(new SingleHit(1, 2, 2.0F, true, 20));
        hitsTen.add(new SingleHit(10, 20, 4.0F, true, 10));
        hitsTen.add(new SingleHit(10, 10, 4.0F, true, 10));
        hitsTen.add(new SingleHit(10, 40, 4.0F, false, 10));

        c.storeSingle(name,hitsOne);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertFalse(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(), 1);

        c.storeSingle(name,hitsTen);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertTrue(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(),2);

        List<SingleHit> hits = new ArrayList<SingleHit>();
        hits.addAll(hitsOne);
        hits.addAll(hitsTen);
        Collections.sort(hits);
        
        assertEquals(name + " chrom 1 has 4 hits (either strand)", c.getCount(name, 1, false,null,null,null,null,null), 4);
        assertEquals(name + " chrom 10 has 4 hits (either strand)", c.getCount(name, 10, false,null,null,null,null,null),3);
        assertEquals(name + " chrom 1 has 4 hits (+ strand)", c.getCount(name, 1, false,null,null,null,null,true),4);
        assertEquals(name + " chrom 10 has 4 hits (+ strand)", c.getCount(name, 10, false,null,null,null,null,true),2);
        assertEquals(name + " chrom 1 has 0 hits (- strand)", c.getCount(name, 1, false,null,null,null,null,false), 0);
        assertEquals(name + " chrom 10 has 1 hits (- strand)", c.getCount(name, 10, false,null,null,null,null,false), 1);

        List<SingleHit> retrieved = c.getSingleHits(name,1,null,null,null,null);
        assertEquals(name + " get chr 1 hits", retrieved.size(),4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits pos " + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,true);
        assertEquals(name + " get chr 1 hits +", retrieved.size(),4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits + pos " + i  + i + ": " + retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,false);
        assertEquals(name + " get chr 1 hits -", retrieved.size(),0);

        retrieved = c.getSingleHits(name,10,null,null,null,null);
        assertEquals(name + " get chr 10 hits", retrieved.size(),3);
        for (int i = 0; i < 3; i++) {
            assertTrue(name + " get chr 10 hits pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,true);
        assertEquals(name + " get chr 10 hits +", retrieved.size(),2);
        for (int i = 0; i < 2; i++) {
            assertTrue(name + " get chr 10 hits + pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,false);
        assertEquals(name + " get chr 10 hits -", retrieved.size(),1);
        for (int i = 0; i < 1; i++) {
            assertTrue(name + " get chr 10 hits - pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+6)));
        }

        c.close();
    }
    @Test public void storeUnsortedSingle() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeUnsortedSingle";

        ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
        hits.add(new SingleHit(10, 10, 4.0F, true, 10));
        hits.add(new SingleHit(1, 2, 2.0F, true, 20));
        hits.add(new SingleHit(10, 40, 4.0F, false, 10));
        hits.add(new SingleHit(10, 20, 4.0F, true, 10));
        hits.add(new SingleHit(1, 1000, 4.0F, true, 10));
        hits.add(new SingleHit(1, 3, 3.0F, true, 30));
        hits.add(new SingleHit(1, 1, 1.0F, true, 10));

        c.storeSingle(name,hits);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertTrue(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(), 2);
        
        Collections.sort(hits);

        assertEquals(name + " chrom 1 has 4 hits (either strand)", c.getCount(name, 1, false,null,null,null,null,null), 4);
        assertEquals(name + " chrom 10 has 4 hits (either strand)", c.getCount(name, 10, false,null,null,null,null,null), 3);
        assertEquals(name + " chrom 1 has 4 hits (+ strand)", c.getCount(name, 1, false,null,null,null,null,true), 4);
        assertEquals(name + " chrom 10 has 4 hits (+ strand)", c.getCount(name, 10, false,null,null,null,null,true), 2);
        assertEquals(name + " chrom 1 has 0 hits (- strand)", c.getCount(name, 1, false,null,null,null,null,false), 0);
        assertEquals(name + " chrom 10 has 1 hits (- strand)", c.getCount(name, 10, false,null,null,null,null,false), 1);

        List<SingleHit> retrieved = c.getSingleHits(name,1,null,null,null,null);
        assertEquals(name + " get chr 1 hits", retrieved.size(), 4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,true);
        assertEquals(name + " get chr 1 hits +", retrieved.size(), 4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 1 hits + pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,false);
        assertEquals(name + " get chr 1 hits -", retrieved.size(), 0);

        retrieved = c.getSingleHits(name,10,null,null,null,null);
        assertEquals(name + " get chr 10 hits", retrieved.size(), 3);
        for (int i = 0; i < 3; i++) {
            assertTrue(name + " get chr 10 hits pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,true);
        assertEquals(name + " get chr 10 hits +", retrieved.size(), 2);
        for (int i = 0; i < 2; i++) {
            assertTrue(name + " get chr 10 hits + pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+4)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,false);
        assertEquals(name + " get chr 10 hits -", retrieved.size(), 1);
        for (int i = 0; i < 1; i++) {
            assertTrue(name + " get chr 10 hits - pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+6)));
        }

        c.close();
    }
    @Test public void storeDouble() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "storeDouble";

        ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
        hits.add(new SingleHit(10, 10, 4.0F, true, 10));
        hits.add(new SingleHit(10, 40, 4.0F, false, 10));
        hits.add(new SingleHit(10, 20, 4.0F, true, 10));
        hits.add(new SingleHit(1, 1000, 4.0F, true, 10));
        hits.add(new SingleHit(1, 3, 3.0F, true, 30));
        hits.add(new SingleHit(1, 1, 1.0F, true, 10));
        hits.add(new SingleHit(1, 2, 2.0F, true, 20));

        c.storeSingle(name,hits);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertTrue(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(), 2);

        assertEquals(name + " chrom 1 has 8 hits (either strand)", c.getCount(name, 1, false,null,null,null,null,null), 4);
        assertEquals(name + " chrom 10 has 6 hits (either strand)", c.getCount(name, 10, false,null,null,null,null,null), 3);
        assertEquals(name + " chrom 1 has 8 hits (+ strand)", c.getCount(name, 1, false,null,null,null,null,true), 4);
        assertEquals(name + " chrom 10 has 4 hits (+ strand)", c.getCount(name, 10, false,null,null,null,null,true), 2);
        assertEquals(name + " chrom 1 has 0 hits (- strand)", c.getCount(name, 1, false,null,null,null,null,false), 0);
        assertEquals(name + " chrom 10 has 2 hits (- strand)", c.getCount(name, 10, false,null,null,null,null,false), 1);

        c.storeSingle(name,hits);        
        assertTrue(name + " exists",c.exists(name));
        assertTrue(name + " has chrom 1",c.getChroms(name, false,false).contains(1));
        assertTrue(name + " has chrom 10",c.getChroms(name, false,false).contains(10));
        assertEquals(name + " has two chroms",c.getChroms(name,false,false).size(), 2);

        hits.addAll(hits);
        Collections.sort(hits);
        
        assertEquals(name + " chrom 1 has 8 hits (either strand)", c.getCount(name, 1, false,null,null,null,null,null), 8);
        assertEquals(name + " chrom 10 has 6 hits (either strand)", c.getCount(name, 10, false,null,null,null,null,null), 6);
        assertEquals(name + " chrom 1 has 8 hits (+ strand)", c.getCount(name, 1, false,null,null,null,null,true), 8);
        assertEquals(name + " chrom 10 has 4 hits (+ strand)", c.getCount(name, 10, false,null,null,null,null,true), 4);
        assertEquals(name + " chrom 1 has 0 hits (- strand)", c.getCount(name, 1, false,null,null,null,null,false), 0);
        assertEquals(name + " chrom 10 has 2 hits (- strand)", c.getCount(name, 10, false,null,null,null,null,false), 2);

        List<SingleHit> retrieved = c.getSingleHits(name,1,null,null,null,null);
        assertEquals(name + " get chr 1 hits", retrieved.size(), 8);
        for (int i = 0; i < 8; i++) {
            assertTrue(name + " get chr 1 hits pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,true);
        assertEquals(name + " get chr 1 hits +", retrieved.size(), 8);
        for (int i = 0; i < 8; i++) {
            assertTrue(name + " get chr 1 hits + pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i)));
        }
        retrieved = c.getSingleHits(name,1,null,null,null,false);
        assertEquals(name + " get chr 1 hits -", retrieved.size(), 0);

        retrieved = c.getSingleHits(name,10,null,null,null,null);
        assertEquals(name + " get chr 10 hits", retrieved.size(), 6);
        for (int i = 0; i < 6; i++) {
            assertTrue(name + " get chr 10 hits pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+8)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,true);
        assertEquals(name + " get chr 10 hits +", retrieved.size(), 4);
        for (int i = 0; i < 4; i++) {
            assertTrue(name + " get chr 10 hits + pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+8)));
        }
        retrieved = c.getSingleHits(name,10,null,null,null,false);
        assertEquals(name + " get chr 10 hits -", retrieved.size(), 2);
        for (int i = 0; i < 1; i++) {
            assertTrue(name + " get chr 10 hits - pos " + i+ retrieved.get(i) + " vs " + hits.get(i), 
                       retrieved.get(i).equals(hits.get(i+12)));
        }

        c.close();
    }    
    @Test public void testACL() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testACL";
        
        int chrom = 60;
        ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
        hits.add(new SingleHit(chrom, 1, 1.0F, true, 10));
        c.storeSingle(name,hits);

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
            other.getChroms(name,false,false);
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
        other.getChroms(name,false,false);

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
        int chrom = 10101010;
        int n = 10000;
        List<SingleHit> hits = new ArrayList<SingleHit>();
        for (int i = 0; i < n; i++) {
            hits.add(new SingleHit(chrom, (int)Math.round(Math.random() * Integer.MAX_VALUE), 1.0F, true, 20));
        }
        c.storeSingle(name,hits);
        
        Collections.sort(hits);
        for (int q = 0; q < 30; q++) {
            int start = (int)Math.round(Math.random() * Integer.MAX_VALUE);
            int end = start + (int)(Math.round(Math.random() * Integer.MAX_VALUE) % (Integer.MAX_VALUE - start));
            int count = 0;
            for (int i = 0; i < hits.size(); i++) {
                if (hits.get(i).pos >= start && hits.get(i).pos <= end) {
                    count++;
                }
            }
            assertTrue(String.format("%s %d != %d", name, count, c.getCount(name,chrom,false,start,end,null,null,null)), count == c.getCount(name,chrom,false,start,end,null,null,null));
            int back[] = c.getPositions(name,chrom,false,start,end,null,null,null);
            int j = 0;
            int k = 0;
            while (j < hits.size() && hits.get(j).pos < start) { j++;}
            while (j < hits.size() && hits.get(j).pos <= end) {
                assertTrue(hits.get(j++).pos == back[k++]);
            }
        }
        assertTrue(c.getCount(name,chrom,false,0, Integer.MAX_VALUE - 1,null,null,null) == n);
        System.err.println("Done with testRangeQuery");
        c.close();
    }
    @Test public void testHistogram() throws IOException, ClientException {
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testHistogram";
        int chrom = 60;
        ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
        hits.add(new SingleHit(chrom, 1, 10.0F, true, 10));
        hits.add(new SingleHit(chrom, 1, 10.0F, false, 10));
        hits.add(new SingleHit(chrom, 1, 10.0F, true, 10));
        hits.add(new SingleHit(chrom, 10, 10.0F, false, 10));
        hits.add(new SingleHit(chrom, 10, 20.0F, true, 20));
        hits.add(new SingleHit(chrom, 11, 20.0F, false, 20));
        hits.add(new SingleHit(chrom, 12, 20.0F, true, 20));
        hits.add(new SingleHit(chrom, 25, 20.0F, false, 20));
        
        c.storeSingle(name,hits);
        
        TreeMap<Integer,Integer> map = c.getHistogram(name,chrom,false,false,10,0,100,null,null);
        System.out.println(name + " MAP1 IS " + map);
        assertEquals(name + " basic size",map.size(),3);
        assertEquals(name + " basic 5",(int)map.get(5),3);
        assertEquals(name + " basic 15",(int)map.get(15),4);
        assertEquals(name + " basic 25",(int)map.get(25),1);

        map = c.getHistogram(name,chrom,false,false,10,0,100,15F,null);
        System.out.println(name + " MAP2 IS " + map);
        assertEquals(name + " minweight size",map.size(),2);
        assertEquals(name + " minweight 15",(int)map.get(15),3);
        assertEquals(name + " minweight 25",(int)map.get(25),1);

        map = c.getHistogram(name,chrom,false,false,10,0,100,null,true);
        System.out.println(name + " MAP3 IS " + map);
        assertEquals(name + " basic size",map.size(),2);
        assertEquals(name + " basic 5",(int)map.get(5),2);
        assertEquals(name + " basic 15",(int)map.get(15),2);
        c.close();
     }
    @Test public void testHistogram2() throws IOException, ClientException {
        int MAXVALUE = 10000;
        float MAXWEIGHT = 4f;
        Client c = new Client(hostname, portnum, user, passwd);
        String name = "testHistogram2";
        for (int q = 0; q < 10; q++) {
            int chrom = q;
            ArrayList<SingleHit> hits = new ArrayList<SingleHit>();
            for (int i = 0; i < MAXVALUE; i++) {
                hits.add(new SingleHit(chrom, (int)Math.round(Math.random() * MAXVALUE), (float)(Math.random() * MAXWEIGHT), false, 10));
            }
            c.storeSingle(name,hits);           
            Collections.sort(hits);

            for (int q2 = 0; q2 < 10; q2++) {
                int start = (int)Math.round(Math.random() * MAXVALUE);
                int end = start + (int)(Math.round(Math.random() * MAXVALUE) % (MAXVALUE - start));
                float weight = (float)Math.random() * MAXWEIGHT;
                int binsize = 10 + (int)Math.round(Math.random() * 40);
                if (end < start) {
                    throw new RuntimeException("end < start");
                }
                Map<Integer,Integer> plainhist = c.getHistogram(name,chrom,false,false,binsize,start,end,null,null);
                Map<Integer,Integer> minweighthist = c.getHistogram(name,chrom,false,false,binsize,start,end,weight,null);
                Map<Integer,Float> weighthist = c.getWeightHistogram(name,chrom,false,false,binsize,start,end,null,null);
                Map<Integer,Float> minwwhist = c.getWeightHistogram(name,chrom,false,false,binsize,start,end,weight,null);

                Map<Integer,Integer> myplain = new HashMap<Integer,Integer>();
                Map<Integer,Integer> myminw = new HashMap<Integer,Integer>();
                Map<Integer,Float> myw = new HashMap<Integer,Float>();
                Map<Integer,Float> myww = new HashMap<Integer,Float>();
                for (int i = 0; i < hits.size(); i++) {
                    int bin = ((hits.get(i).pos - start) / binsize) * binsize + binsize / 2 + start;
                    if (hits.get(i).pos >= start && hits.get(i).pos <= end) {
                        if (myplain.containsKey(bin)) {
                            myplain.put(bin,1 + myplain.get(bin));
                        } else {
                            myplain.put(bin,1);
                        }
                        if (myw.containsKey(bin)) {
                            myw.put(bin,hits.get(i).weight + myw.get(bin));
                        } else {
                            myw.put(bin,hits.get(i).weight);
                        }                        
                        if (hits.get(i).weight >= weight) {
                            if (myminw.containsKey(bin)) {
                                myminw.put(bin,1 + myminw.get(bin));
                            } else {
                                myminw.put(bin,1);
                            }
                            if (myww.containsKey(bin)) {
                                myww.put(bin,hits.get(i).weight + myww.get(bin));
                            } else {
                                myww.put(bin,hits.get(i).weight);
                            }                            
                        }
                    }
                }
                System.err.println("MyW " + myw);
                System.err.println("weighthist " +weighthist);
                assertEquals(myplain.size(),plainhist.size());
                assertEquals(myminw.size(),minweighthist.size());
                assertEquals(myw.size(),weighthist.size());
                assertEquals(myww.size(),minwwhist.size());
                for (int bin : myplain.keySet()) {
                    assertEquals(String.format("plain bin %d : %d != %d", bin, myplain.get(bin), plainhist.get(bin)), 
                               myplain.get(bin),plainhist.get(bin));
                }
                for (int bin : myminw.keySet()) {
                    assertEquals(String.format("minw bin %d : %d != %d", bin, myminw.get(bin), minweighthist.get(bin)), 
                               myminw.get(bin),minweighthist.get(bin));                   
                }
                for (int bin : myw.keySet()) {
                    assertEquals(String.format("weight bin %d : %f != %f", bin, myw.get(bin), weighthist.get(bin)), 
                                 myw.get(bin),weighthist.get(bin), .0001);
                }
                for (int bin : myww.keySet()) {
                    assertEquals(String.format("ww bin %d : %f != %f", bin, myww.get(bin), minwwhist.get(bin)), 
                                 myww.get(bin) , minwwhist.get(bin), .0001);
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

        org.junit.runner.JUnitCore.main("edu.mit.csail.cgs.projects.readdb.TestReadDB");
    }
}