package edu.mit.csail.cgs.projects.readdb;

import java.util.TreeMap;
import java.io.IOException;

public interface ReadOnlyClient {

    public TreeMap<Integer,Integer> getHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException;
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, int chromid, boolean paired, boolean doReadExtension, int binsize, Integer start, Integer stop, Float minWeight, Boolean plusStrand) throws IOException, ClientException;
    public void close();
}