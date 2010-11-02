package edu.mit.csail.cgs.projects.readdb;

import java.io.IOException;
import java.util.TreeMap;
import java.util.Map;
import java.util.ArrayList;
import java.util.Collection;

public class Aggregator implements ReadOnlyClient {

    private ArrayList<ReadOnlyClient> clients;
    public Aggregator() {
        clients = new ArrayList<ReadOnlyClient>();        
    }
    public Aggregator(Collection<? extends ReadOnlyClient> c) {
        clients = new ArrayList<ReadOnlyClient>();        
        clients.addAll(c);
    }

    /**
     * merges a into b, unless b is null in which case returns a
     */
    public static TreeMap<Integer,Integer> mergeHistogramsII(TreeMap<Integer,Integer> a, TreeMap<Integer,Integer> b) {
        if (a == null) {
            return b;
        } else if (b == null) {
            return a;
        }
        for (int k : a.keySet()) {
            if (b.containsKey(k)) {
                b.put(k, b.get(k) + a.get(k));
            } else {
                b.put(k,a.get(k));
            }
        }
        return b;
    }
    public static TreeMap<Integer,Float> mergeHistogramsFF(TreeMap<Integer,Float> a, TreeMap<Integer,Float> b) {
        if (a == null) {
            return b;
        } else if (b == null) {
            return a;
        }
        for (int k : a.keySet()) {
            if (b.containsKey(k)) {
                b.put(k, b.get(k) + a.get(k));
            } else {
                b.put(k,a.get(k));
            }
        }
        return b;
    }
    public static TreeMap<Integer,Float> mergeHistogramsIF(TreeMap<Integer,Integer> a, TreeMap<Integer,Float> b) {
        if (a == null) {
            return b;
        } else if (b == null) {
            b = new TreeMap<Integer,Float>();
        }
        for (int k : a.keySet()) {
            if (b.containsKey(k)) {
                b.put(k, b.get(k) + a.get(k));
            } else {
                b.put(k,(float)a.get(k));
            }
        }
        return b;
    }
    

    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        if (clients.size() == 0) {
            return new TreeMap<Integer,Integer>();
        }
        TreeMap<Integer,Integer> output = clients.get(0).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
        for (int i = 1; i < clients.size(); i++) {
            TreeMap<Integer,Integer> o = clients.get(i).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
            output = mergeHistogramsII(o,output);
        }
        return output;
    }


    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        if (clients.size() == 0) {
            return new TreeMap<Integer,Float>();
        }
        TreeMap<Integer,Float> output = clients.get(0).getWeightHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
        for (int i = 1; i < clients.size(); i++) {
            TreeMap<Integer,Float> o = clients.get(i).getWeightHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
            output = mergeHistogramsFF(o,output);
        }
        return output;
    }

    public void close() {
        for (ReadOnlyClient c : clients) {
            c.close();
        }

    }
}