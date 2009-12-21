/*
 * Created on Oct 6, 2006
 */
package edu.mit.csail.cgs.utils.parsing.ncbi;

import java.util.*;
import java.io.*;

/**
 * @author tdanford
 */
public class GeneInfoParser implements Iterator {
    
    public static void printNameMap(File f, String taxid) throws IOException { 
        PrintStream ps = new PrintStream(new FileOutputStream(new File(taxid + "_names.txt")));
        GeneInfoParser parser = new GeneInfoParser(f);
        while(parser.hasNext()) { 
            Entry e = (Entry)parser.next();
            if(e.tax_id.equals(taxid)) { 
                ps.print(e.geneID +"\t");
				ps.print(e.symbol);
                for(String syn : e.synonyms) { 
					ps.print("|" + syn);
                }
                ps.println();
            }
        }
        ps.close();
    }
    
    public static void main(String[] args) { 
        File infile = new File(args[0]);
        for(int i = 1; i < args.length; i++) { 
            try {
                printNameMap(infile, args[i]);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
    
    private String nextLine;
    private BufferedReader br;
    
    public GeneInfoParser(File f) throws IOException {
        br = new BufferedReader(new FileReader(f));
        br.readLine();
        nextLine = br.readLine();
    }
    
    public boolean hasNext() { 
        return nextLine != null;
    }
    
    public Object next() { 
        String line = nextLine;
        try { 
            nextLine = br.readLine();
            if(nextLine == null) { br.close(); }
        } catch(IOException ie) { 
            ie.printStackTrace();
            nextLine = null;
        }
        return new Entry(line);
    }
    
    public static class Entry { 
        public String tax_id, geneID, symbol, locusTag;
        public LinkedList<String> synonyms;
        public String dbxref, chromosome;

        public Entry(String line) { 
            String[] array = line.split("\t");
            tax_id = array[0];
            geneID = array[1];
            symbol = array[2];
            locusTag = array[3];
            synonyms = new LinkedList<String>();
			if(!array[4].equals("-")) { 
				String[] synarray = array[4].split("\\|");
				for(int i = 0; i < synarray.length; i++) { synonyms.add(synarray[i]); }
			}
            dbxref = array[5];
            chromosome = array[6];
        }
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
