package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.mit.csail.cgs.utils.Closeable;

public class HMMReads implements Closeable {

    private Client client;
    private double subsample;

    public HMMReads() throws ClientException, IOException {
        client = new Client();
        subsample = 2;
    }
    public void subSample(double d) {
        subsample = d;
    }
    public double subSample() {return subsample;}
    public ReadCounts getReadCounts(Region region,
                                    Collection<ChipSeqAlignment> alignments) throws IOException, ClientException {
        int output[] = new int[region.getWidth()];
        int regionstart = region.getStart();
        for (int i = 0; i < output.length; i++) {
            output[i] = 0;
        }
        for (ChipSeqAlignment a : alignments) {
            TreeMap<Integer,Integer> m = client.getHistogram(Integer.toString(a.getDBID()),
                                                             region.getGenome().getChromID(region.getChrom()),
                                                             false,
                                                             false,
                                                             1, // binsize
                                                             10, // dedup
                                                             region.getStart(),
                                                             region.getEnd(),
                                                             null,
                                                             true,
                                                             true);
            for (int pos : m.keySet()) {
                output[pos - regionstart] += m.get(pos);
            }
        }
        if (subsample < 1) {
            Binomial binomial = new Binomial(5, .5, new DRand());
            int newout[] = new int[output.length];
            for (int i = 0; i < newout.length; i++) {
                if (output[i] > 0) {
                    newout[i] = binomial.nextInt(output[i], subsample);
                } else {
                    newout[i] = 0;
                }



            }
            output = newout;
        }

        return new ReadCounts(output,region);
    }

    public void close() {
        try {
            client.close();
            client = null;
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    public boolean isClosed() {
        return client == null;
    }


}