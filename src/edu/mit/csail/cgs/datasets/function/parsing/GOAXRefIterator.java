/*
 * Created on Nov 9, 2006
 */
package edu.mit.csail.cgs.datasets.function.parsing;

import java.util.*;
import java.io.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class GOAXRefIterator implements Iterator<GOAXRefLine> {
    
    public static void main(String[] args) { 
        File f = args.length > 0 ? new File(args[0]) : 
            new File("C:\\Documents and Settings\\tdanford\\Desktop\\go_annotations\\human.xrefs");
        try {
            GOAXRefIterator itr = new GOAXRefIterator(f);
            while(itr.hasNext()) { 
                GOAXRefLine line = itr.next();
                System.out.print(line.getDatabaseID() + " --> ");
                for(Pair<String,String> pair : line.getEntrezIDs()) {
                	String entrez = pair.getFirst();
                	String symbol = pair.getLast();
                    System.out.print(entrez + " ");
                }
                System.out.println();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private GOAXRefLine nextLine;
    private BufferedReader br;
    
    public GOAXRefIterator(File f) throws IOException { 
        br = new BufferedReader(new FileReader(f));
        getNextLine();
    }
    
    private void getNextLine() { 
        String line = null;
        try { 
            do { 
                line = br.readLine();
                if(line != null) { line = line.trim(); }
            } while(line != null && (line.startsWith("#") || line.startsWith("!")));
        } catch(IOException ie) { 
            line = null;
            ie.printStackTrace(System.err);
        }
        
        if(line != null) {
            nextLine = new GOAXRefLine(line);
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

    public GOAXRefLine next() {
        GOAXRefLine ret = nextLine;
        getNextLine();
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
