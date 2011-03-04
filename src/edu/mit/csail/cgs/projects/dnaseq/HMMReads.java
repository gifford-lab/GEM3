package edu.mit.csail.cgs.projects.dnaseq;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.utils.Closeable;

public class HMMReads implements Closeable {

    private Client client;

    public HMMReads() throws ClientException, IOException {
        client = new Client();
    }
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