package edu.mit.csail.cgs.projects.readdb;

import java.util.TreeMap;
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

    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        if (clients.size() == 0) {
            return new TreeMap<Integer,Integer>();
        }
        TreeMap<Integer,Integer> output = clients.get(0).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
        for (int i = 1; i < clients.size(); i++) {
            TreeMap<Integer,Integer> o = clients.get(i).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
            for (int k : o.keySet()) {
                if (output.containsKey(k)) {
                    output.put(o, output.get(k) + o.get(k));
                } else {
                    output.put(o.get(k));
                }
            }
        }
        return output;
    }


    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException {
        if (clients.size() == 0) {
            return new TreeMap<Integer,Integer>();
        }
        TreeMap<Integer,Float> output = clients.get(0).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
        for (int i = 1; i < clients.size(); i++) {
            TreeMap<Integer,Float> o = clients.get(i).getHistogram(alignid,chromid,paired,doReadExtension,binsize,start,stop,minWeight,plusStrand);
            for (int k : o.keySet()) {
                if (output.containsKey(k)) {
                    output.put(o, output.get(k) + o.get(k));
                } else {
                    output.put(o.get(k));
                }
            }
        }
        return output;
    }

    public void close() {
        for (ReadOnlyClient c : clients) {
            c.close();
        }

    }
}