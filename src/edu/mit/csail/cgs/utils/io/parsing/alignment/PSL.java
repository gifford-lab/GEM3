package edu.mit.csail.cgs.utils.io.parsing.alignment;

import java.io.*;
import java.util.Iterator;

public class PSL implements Iterator<PSLHit> {

    public static final int match = 0,
        mismatch = 1,
        repmatch = 2,
        n = 3,
        qgapcount = 4,
        qgapbases = 5,
        tgapcount = 6,
        tgapbases = 7,
        strand = 8,
        qname = 9,
        qsize = 10,
        qstart = 11,
        qend = 12,
        tname = 13,
        tsize = 14,
        tstart = 15,
        tend = 16,
        blockcount = 17,
        blocksizes = 18,
        qstarts = 19,
        tstarts = 20;
    
    private BufferedReader reader;
    private String nextLine = null;
    
    public PSL(File f) throws IOException { 
    	reader = new BufferedReader(new FileReader(f));
        nextLine = reader.readLine();
        while (nextLine != null &&
               !nextLine.matches("^\\d.*")) {
            nextLine = reader.readLine();
        }
    }
    
    public PSL (BufferedReader reader) throws IOException {
        this.reader = reader;
        nextLine = reader.readLine();
        while (nextLine != null &&
               !nextLine.matches("^\\d.*")) {
            nextLine = reader.readLine();
        }

    }
    public void remove() {throw new UnsupportedOperationException("Can't remove from PSL");}
    public boolean hasNext() {return nextLine != null;}
    public PSLHit next() {
        nextLine = nextLine.trim();
        String pieces[] = nextLine.split("\\t");
        try {
            nextLine = reader.readLine();
        } catch (IOException e) {
            nextLine = null;
            throw new RuntimeException("can't read" + e.toString(),e);
        }

        PSLHit hit = new PSLHit();
        hit.match = Integer.parseInt(pieces[match]);
        hit.mismatch = Integer.parseInt(pieces[mismatch]);
        hit.repmatch = Integer.parseInt(pieces[repmatch]);
        hit.n = Integer.parseInt(pieces[n]);
        hit.qgapcount = Integer.parseInt(pieces[qgapcount]);
        hit.tgapcount = Integer.parseInt(pieces[tgapcount]);
        hit.qgapbases = Integer.parseInt(pieces[qgapbases]);
        hit.tgapbases = Integer.parseInt(pieces[tgapbases]);
        hit.strand = pieces[strand].charAt(0);
        hit.qname = pieces[qname];
        hit.tname = pieces[tname];
        hit.qsize = Integer.parseInt(pieces[qsize]);
        hit.qstart = Integer.parseInt(pieces[qstart]);
        hit.qend = Integer.parseInt(pieces[qend]);
        hit.tsize = Integer.parseInt(pieces[tsize]);
        hit.tstart = Integer.parseInt(pieces[tstart]);
        hit.tend = Integer.parseInt(pieces[tend]);
        int blocks = Integer.parseInt(pieces[blockcount]);
        hit.blockSizes = new int[blocks];
        hit.qStarts = new int[blocks];
        hit.tStarts = new int[blocks];
        String sub[] = pieces[blocksizes].split(",");
        for (int i = 0; i < blocks; i++) {
            hit.blockSizes[i] = Integer.parseInt(sub[i]);
        }
        sub = pieces[qstarts].split(",");
        for (int i = 0; i < blocks; i++) {
            hit.qStarts[i] = Integer.parseInt(sub[i]);
        }
        sub = pieces[tstarts].split(",");
        for (int i = 0; i < blocks; i++) {
            hit.tStarts[i] = Integer.parseInt(sub[i]);
        }
        return hit;
    }


}
