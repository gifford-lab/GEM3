package edu.mit.csail.cgs.projects.readdb;

import java.nio.*;

public class ReIndexHits {
    /**
     * Regenerates the header files with a new pagesize
     */
    public static void main(String args[]) throws Exception {
        int pagesize = Integer.parseInt(args[0]);
        for (int i = 1; i < args.length; i++) {
            String hitsFname = args[i].replaceAll("index$", "hits");
            String weightsFname = args[i].replaceAll("index$", "weights");
            Hits hits = new Hits(hitsFname, weightsFname);
            IntBuffer h = hits.getPositionsBuffer().ib;
            Header newheader = new Header(h, pagesize);
            newheader.writeIndexFile(args[i]);
        }

    }

}