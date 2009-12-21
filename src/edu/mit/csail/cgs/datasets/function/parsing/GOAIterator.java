/*
 * Created on Nov 9, 2006
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class GOAIterator implements Iterator<GOALine> {
    
    public static Map<String,Set<String>> createMap(File f, int i1, int i2) {  
        Map<String,Set<String>> map = new HashMap<String,Set<String>>();
        try {
            GOAIterator parser = new GOAIterator(f);
            while(parser.hasNext()) { 
                GOALine line = parser.next();
                String k = line.array[i1];
                String v = line.array[i2];
                if(!map.containsKey(k)) { map.put(k, new HashSet<String>()); }
                map.get(k).add(v);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return map;
    }
    
    private GOALine nextLine;
    private BufferedReader br;
    
    public GOAIterator(File f) throws IOException { 
        br = new BufferedReader(new FileReader(f));
        getNextLine();
    }
    
    private void getNextLine() { 
        String line = null;
        try { 
            do { 
                line = br.readLine();
            } while(line != null && line.startsWith("!"));
        } catch(IOException ie) { 
            line = null;
            ie.printStackTrace(System.err);
        }
        
        if(line != null) { 
            nextLine = new GOALine(line);
        } else { 
            nextLine = null;
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public boolean hasNext() {
        return nextLine != null;
    }

    public GOALine next() {
        GOALine ret = nextLine;
        getNextLine();
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
