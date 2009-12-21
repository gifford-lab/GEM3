package edu.mit.csail.cgs.projects.readdb;

import java.util.*;
import java.io.IOException;

public interface ReadOnlyClient {

    public boolean exists(String alignid) throws IOException;
    public Set<String> getChroms(String alignid) throws IOException, ClientException;
    public int getCount(String alignid) throws IOException, ClientException;
    public double getWeight(String alignid) throws IOException, ClientException;
    public int getCount(String alignid, String chromid)  throws IOException, ClientException;
    public double getWeight(String alignid, String chromid) throws IOException, ClientException;
    public int getCountRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException;
    public int getCountRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException;
    public double getWeightRange(String alignid, String chromid, int start, int stop)  throws IOException, ClientException;
    public double getWeightRange(String alignid, String chromid, int start, int stop, float minweight)  throws IOException, ClientException;
    public int[] getHitsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException;
    public int[] getHitsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException;
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop) throws IOException, ClientException;
    public float[] getWeightsRange(String alignid, String chromid, int start, int stop, float minweight) throws IOException, ClientException;
    public int[] getHits(String alignid, String chromid)  throws IOException, ClientException;
    public float[] getWeights(String alignid, String chromid)  throws IOException, ClientException;
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException;
    public TreeMap<Integer,Integer> getHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension) throws IOException, ClientException;
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize) throws IOException, ClientException;
    public TreeMap<Integer,Float> getWeightHistogram(String alignid, String chromid, int start, int stop, int binsize, float minweight, int readExtension)  throws IOException, ClientException;
    public void close();
}