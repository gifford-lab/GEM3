package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.IOException;

/**
 * Aggregator performs operations on a set of alignments
 * and returns results that are aggregated across all of them.
 * You might use this, for example, to retrieve results
 * from all replicates of an experiment.
 */

public class Aggregator implements ReadOnlyClient {

    private ReadOnlyClient client;
    
    public Aggregator(ReadOnlyClient c) {
        client = c;
    }
    /* first section is the ReadOnlyClient methods.
       
       Aggregator is here to
       - keep Client as simple as possible (we could move the Set<Alignid> methods there)
       - allow code to just use Aggregator if it wants to either work on a set of
         alignments or a single alignment
       - it ignores the read-write methods because the semantics of those aren't
         clear across multiple alignments

     */
    public boolean exists(String alignid) throws IOException {return client.exists(alignid);}
    public Set<String> getChroms(String alignid) throws IOException, ClientException {return client.getChroms(alignid);}
    public int getCount(String alignid) throws IOException, ClientException {return client.getCount(alignid);}
    public double getWeight(String alignid) throws IOException, ClientException {return client.getWeight(alignid);}
    public int getCount(String alignid, String chromid)  throws IOException, ClientException {return client.getCount(alignid, chromid);}
    public double getWeight(String alignid, String chromid) throws IOException, ClientException {return client.getWeight(alignid, chromid);}
    public int getCountRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException {return getCountRange(alignid, chromid, start, stop);}
    public int getCountRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {return client.getCountRange(alignid, chromid, start, stop, minweight);}
    public double getWeightRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException {return client.getWeightRange(alignid, chromid, start, stop);}
    public double getWeightRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {return client.getWeightRange(alignid, chromid, start, stop, minweight);}
    public int[] getHitsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException {return client.getHitsRange(alignid, chromid, start, stop);}
    public int[] getHitsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {return client.getHitsRange(alignid, chromid, start, stop,minweight);}
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException {return client.getWeightsRange(alignid, chromid, start, stop);}
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {return client.getWeightsRange(alignid, chromid, start, stop, minweight);}
    public int[] getHits(String alignid, String chromid)  throws IOException, ClientException {return client.getHits(alignid, chromid);}
    public float[] getWeights(String alignid, String chromid)  throws IOException, ClientException {return client.getWeights(alignid, chromid);}
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {return client.getHistogram(alignid, chromid, start, stop, binsize);}
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) throws IOException, ClientException {return client.getHistogram(alignid, chromid, start, stop, binsize, minweight, readExtension);}
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {return client.getWeightHistogram(alignid, chromid, start, stop, binsize);}
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension)  throws IOException, ClientException {return client.getWeightHistogram(alignid, chromid, start, stop, binsize, minweight, readExtension);}
    // second part is the aggregator methods
    public Set<String> getChroms(Set<String> alignid) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getChroms(getOnly(alignid));
        }
        Set<String> output = new HashSet<String>();
        for (String s : alignid) {
            output.addAll(client.getChroms(s));
        }
        return output;
    }
    public int getCount(Set<String> alignid) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getCount(getOnly(alignid));
        }
        int c = 0;
        for (String s : alignid) {
            c += client.getCount(s);
        }
        return c;
    }
    public double getWeight(Set<String> alignid) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeight(getOnly(alignid));
        }
        double d = 0;
        for (String s : alignid) {
            d += client.getWeight(s);
        }
        return d;
    }
    public int getCount(Set<String> alignid, String chromid)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getCount(getOnly(alignid), chromid);
        }
        int c = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                c += client.getCount(s,chromid);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return c;
    }
    public int getCountRange(Set<String> alignid, String chromid, int start, int stop)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getCountRange(getOnly(alignid), chromid, start, stop);
        }
        int c = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                c += client.getCountRange(s,chromid,start,stop);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return c;
    }
    public int getCountRange(Set<String> alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getCountRange(getOnly(alignid), chromid, start, stop, minweight);
        }
        int c = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                c += client.getCountRange(s,chromid,start,stop,minweight);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return c;

    }
    public double getWeight(Set<String> alignid, String chromid) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeight(getOnly(alignid), chromid);
        }
        double d = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                d += client.getWeight(s,chromid);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return d;
    }
    public double getWeightRange(Set<String> alignid, String chromid, int start, int stop)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightRange(getOnly(alignid), chromid, start, stop);
        }
        double d = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                d += client.getWeightRange(s,chromid, start, stop);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return d;
    }
    public double getWeightRange(Set<String> alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightRange(getOnly(alignid), chromid, start, stop, minweight);
        }
        double d = 0;
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                d += client.getWeightRange(s,chromid, start, stop, minweight);
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return d;
    }
    private int[] mergeHits(int[][] positions, int count) {
        int[] output = new int[count];
        int[] pointers = new int[positions.length];
        for (int i = 0; i < pointers.length; i++) {
            pointers[i] = 0;
        }
        boolean done = false;
        int index = 0;
        while (!done) {
            int minj = -1, minval = Integer.MAX_VALUE;
            int j = 0;
            while (j < pointers.length) {
                if (pointers[j] < positions[j].length && positions[j][pointers[j]] < minval) {
                    minval = positions[j][pointers[j]];
                    minj = j;
                }
                j++;
            }
            if (minj == -1) {
                done = true;
            } else {
                output[index++] = minval;
                pointers[minj]++;
            }
        }
        return output;
    }
    private float[] mergeHits(int[][] positions, float[][] weights, int count) {
        float[] output = new float[count];
        int[] pointers = new int[positions.length];
        for (int i = 0; i < pointers.length; i++) {
            pointers[i] = 0;
        }
        boolean done = false;
        int index = 0;
        while (!done) {
            int minj = -1, minval = Integer.MAX_VALUE;
            int j = 0;
            while (j < pointers.length) {
                if (pointers[j] < positions[j].length && positions[j][pointers[j]] < minval) {
                    minval = positions[j][pointers[j]];
                    minj = j;
                }
                j++;
            }
            if (minj == -1) {
                done = true;
            } else {
                output[index++] = weights[minj][pointers[minj]];
                pointers[minj]++;
            }
        }
        return output;
    }
    public int[] getHitsRange(Set<String> alignid, String chromid, int start, int stop) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getHitsRange(getOnly(alignid), chromid, start, stop);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j++] = client.getHitsRange(s, chromid, start, stop);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, count);
    }
    public int[] getHitsRange(Set<String> alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getHitsRange(getOnly(alignid), chromid, start, stop, minweight);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j++] = client.getHitsRange(s, chromid, start, stop, minweight);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, count);
    }
    public int[] getHits(Set<String> alignid, String chromid)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getHits(getOnly(alignid), chromid);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j++] = client.getHits(s, chromid);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, count);
    }
    public float[] getWeightsRange(Set<String> alignid, String chromid, int start, int stop) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightsRange(getOnly(alignid), chromid, start, stop);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        float weights[][] = new float[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j] = client.getHitsRange(s, chromid, start,stop);
                weights[j++] = client.getWeightsRange(s, chromid, start, stop);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, weights, count);
    }
    public float[] getWeightsRange(Set<String> alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightsRange(getOnly(alignid), chromid, start, stop, minweight);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        float weights[][] = new float[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j] = client.getHitsRange(s, chromid, start,stop, minweight);
                weights[j++] = client.getWeightsRange(s, chromid, start, stop, minweight);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, weights, count);
    }
    public float[] getWeights(Set<String> alignid, String chromid)  throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeights(getOnly(alignid), chromid);
        }
        int count = 0;
        int positions[][] = new int[alignid.size()][];
        float weights[][] = new float[alignid.size()][];
        boolean anyworked = false;
        int j = 0;
        for (String s : alignid) {
            try {
                positions[j] = client.getHits(s, chromid);
                weights[j++] = client.getWeights(s, chromid);
                count += positions[j-1].length;
                anyworked = true;
            } catch (ClientException e) {
                // ignore it.  not all alignments may have all chromosomes.
            }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return mergeHits(positions, weights, count);
    }
    public TreeMap<Integer,Integer> getHistogram(Set<String> alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getHistogram(getOnly(alignid), chromid, start, stop, binsize);
        }
        TreeMap<Integer, Integer> output = new TreeMap<Integer, Integer>();
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                TreeMap<Integer,Integer> t = client.getHistogram(s, chromid, start, stop, binsize);
                anyworked = true;
                for (int i : t.keySet()) {
                    if (output.containsKey(i)) {
                        output.put(i, output.get(i) + t.get(i));
                    } else {
                        output.put(i, t.get(i));
                    }
                }

            } catch (ClientException e) { }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return output;
    }
    public TreeMap<Integer,Integer> getHistogram(Set<String> alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getHistogram(getOnly(alignid), chromid, start, stop, binsize, minweight, readExtension);
        }
        TreeMap<Integer, Integer> output = new TreeMap<Integer, Integer>();
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                TreeMap<Integer,Integer> t = client.getHistogram(s, chromid, start, stop, binsize, minweight, readExtension);
                anyworked = true;
                for (int i : t.keySet()) {
                    if (output.containsKey(i)) {
                        output.put(i, output.get(i) + t.get(i));
                    } else {
                        output.put(i, t.get(i));
                    }
                }

            } catch (ClientException e) { }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return output;
    }
    public TreeMap<Integer,Float> getWeightHistogram(Set<String> alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightHistogram(getOnly(alignid), chromid, start, stop, binsize);
        }
        TreeMap<Integer, Float> output = new TreeMap<Integer, Float>();
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                TreeMap<Integer,Float> t = client.getWeightHistogram(s, chromid, start, stop, binsize);
                anyworked = true;
                for (int i : t.keySet()) {
                    if (output.containsKey(i)) {
                        output.put(i, output.get(i) + t.get(i));
                    } else {
                        output.put(i, t.get(i));
                    }
                }

            } catch (ClientException e) { }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return output;
    }
    public TreeMap<Integer,Float> getWeightHistogram(Set<String> alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) throws IOException, ClientException {
        if (alignid.size() == 1) {
            return client.getWeightHistogram(getOnly(alignid), chromid, start, stop, binsize, minweight, readExtension);
        }
        TreeMap<Integer, Float> output = new TreeMap<Integer, Float>();
        boolean anyworked = false;
        for (String s : alignid) {
            try {
                TreeMap<Integer,Float> t = client.getWeightHistogram(s, chromid, start, stop, binsize, minweight, readExtension);
                anyworked = true;
                for (int i : t.keySet()) {
                    if (output.containsKey(i)) {
                        output.put(i, output.get(i) + t.get(i));
                    } else {
                        output.put(i, t.get(i));
                    }
                }

            } catch (ClientException e) { }
        }
        if (!anyworked) {
            throw new ClientException("Invalid chromosome " + chromid + " for alignments " + alignid);
        }
        return output;
    }
    public void close() {
        client.close();
    }
    private String getOnly(Set<String> s) {
        String o = null;
        Iterator<String> i = s.iterator();
        o = i.next();
        // shouldn't be any more, but run out the iterator in case
        while (i.hasNext()) {
            i.next();  
        }

        return o;
    }

}