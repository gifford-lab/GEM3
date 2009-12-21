/*
 * Created on Oct 6, 2005
 */
package edu.mit.csail.cgs.utils.parsing;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.Pair;

/**
 * @author tdanford
 *
 * <code>FASTAStream</code> is an iterator over the sequences in a FASTA file.
 * Each sequence is returned as a Pair of Strings: the first String is the
 * name of the sequence and the second String is the sequence itself.
 */
public class FASTAStream implements Iterator<Pair<String,String>>, edu.mit.csail.cgs.utils.Closeable {
    
    private String label;
    private BufferedReader br;
    private String pendingName;

    public FASTAStream(File f) throws IOException {
        br = new BufferedReader(new FileReader(f));
        label = f.getAbsolutePath();
        init();
    }
    public FASTAStream(BufferedReader r) throws IOException {
        br= r;
        label = "inputstream";
        init();
    }
    private void init() throws IOException {
        pendingName = null;
        String line;
        boolean searching = true;
        while(searching && (line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0 && line.charAt(0) == '>') { 
                pendingName = parsePendingName(line);
                searching = false;
            }
        }
    }

    private String parsePendingName(String n) { 
        String base = n.trim();
        if(base.length() > 0 && base.charAt(0) == '>') { 
            base = base.substring(1, base.length()); 
        }
        return base;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.utils.Stream#hasNext()
     */
    public boolean hasNext() {
        return pendingName != null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.utils.Stream#next()
     */
    public Pair<String,String> next() {
        String line = null;
        String name = pendingName;
        pendingName = null;
        StringBuilder sb = new StringBuilder();
        boolean reading = true;
        try { 
            while(reading && (line = br.readLine()) != null) { 
                line = line.trim();
                if(line.length() > 0 && line.charAt(0) == '>') { 
                    pendingName = parsePendingName(line);
                    reading = false;
                } else { 
                    sb.append(line);
                }
            }
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
        
        if(pendingName == null) { close(); }
        Pair<String,String> result = new Pair<String,String>(name, sb.toString());
        sb = null;
        return result;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.hyperdrive.utils.Stream#getDescriptor()
     */
    public String getDescriptor() {
        return "FASTA(" + label + ")";
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        if (isClosed()) {return;}
        try { 
            br.close();
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
        pendingName = null;
        br = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return br == null;
    }
    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a FASTAStream");
    }

}
