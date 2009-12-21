/*
 * Created on Mar 9, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class FileLineExpander implements Expander<File,String> {

    public FileLineExpander() {
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Expander#execute(null)
     */
    public Iterator<String> execute(File f) {
        try {
            return new LazyLineIterator(f);
        } catch (IOException e) {
            e.printStackTrace();
            return new LinkedList<String>().iterator();
        }
    }

    private static class LazyLineIterator implements Iterator<String> {
        
        private BufferedReader br;
        private String nextLine;
        
        public LazyLineIterator(File f) throws IOException { 
            br = new BufferedReader(new FileReader(f));
			if(br == null) { throw new IllegalStateException(); }
            nextLine = br.readLine();
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            return nextLine != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public String next() {
            String ret = nextLine;
            try {
                nextLine = br.readLine();
                if(nextLine == null) { 
                    br.close();
                }
            } catch (IOException e) {
                nextLine = null;
                e.printStackTrace();
            } 
			//System.out.println("\"" + ret + "\"");
            return ret;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        } 
        
    }
}
