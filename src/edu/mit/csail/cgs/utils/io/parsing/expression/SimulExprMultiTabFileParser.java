/*
 * Created on Sep 19, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class SimulExprMultiTabFileParser implements edu.mit.csail.cgs.utils.Closeable {
    
    private String speciesName;
    private Vector<String> exptNames;
    private Vector<LinkedList<ExprProbeValue>> pending;
    private Map<Integer,InnerIterator> iterators;
    private BufferedReader br;

    public SimulExprMultiTabFileParser(File f) throws IOException {
        iterators = new HashMap<Integer,InnerIterator>();
        pending = new Vector<LinkedList<ExprProbeValue>>();
        exptNames = new Vector<String>();
        br = new BufferedReader(new FileReader(f));
        
        String line = br.readLine();
        if(line == null) { throw new IOException("No line/species at start."); }
        int numExpts = 0;
        if(line.indexOf("\t") != -1) { 
            int idx = line.indexOf("\t");
            numExpts = Integer.parseInt(line.substring(0, idx));
            speciesName = line.substring(idx+1, line.length());
        } else { 
            throw new IOException("no <num> <species> line at start.");
        }

        for(int i = 0; i < numExpts; i++) {
            line = br.readLine();
            if(line == null) { throw new IOException("line: " + i); }
            exptNames.add(line);
			System.out.println(i + ": " + exptNames.get(i));
            pending.add(new LinkedList<ExprProbeValue>());
        }
        
        parseNextLine();
    }
    
    public boolean isClosed() { return br == null; }
    public void close() { 
        try { 
            br.close();
            br = null;
        } catch(IOException ie) { 
            ie.printStackTrace(System.err);
        }
    }
    
    public String getSpeciesName() { return speciesName; }
    public int getNumIterators() { return exptNames.size(); }
    public String getExptName(int index) { return exptNames.get(index); }
    
    private void parseNextLine() { 
        if(!isClosed()) { 
            try {
                parseNextLine(br.readLine()); 
            } catch(IOException ie) { 
                ie.printStackTrace(System.err);
            }
        } else { 
            for(int i = 0; i < pending.size(); i++) { 
                if(pending.get(i).isEmpty()) { 
                    pending.set(i, null);
                }
            }            
        }
    }
    
    private void parseNextLine(String line) throws IOException { 
        if(line == null) { 
            if(!isClosed()) { close(); } 
        } else { 
            String[] array = line.split("[\\s]+");
            String name = array[0];
			System.out.print(name + " : "); System.out.flush();
            for(int i = 0; i < exptNames.size(); i++) {
                try { 
                    double val = Double.parseDouble(array[i+1]);
                    ExprProbeValue epv = new ExprProbeValue(exptNames.get(i), name, val);
                    pending.get(i).addLast(epv);
					System.out.print("*"); System.out.flush();
                } catch(NumberFormatException nfe) { 
					System.out.print("-"); System.out.flush();
                }
            }
			System.out.println();
        }
        
    }
    
    public Iterator<ExprProbeValue> getIterator(int index) {
        if(index < 0 || index >= exptNames.size()) { 
            throw new IllegalArgumentException(index + " not in (" + exptNames.size() + ")");
        }
        if(!iterators.containsKey(index)) { iterators.put(index, new InnerIterator(index)); }
        return iterators.get(index);
    }
    
    private class InnerIterator implements Iterator<ExprProbeValue> {
        
        private int index;
        
        public InnerIterator(int index) { this.index = index; }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdrive.utils.Iterator#hasNext()
         */
        public boolean hasNext() {
			while(!isClosed() && pending.get(index).isEmpty()) { 
				parseNextLine(); 
			}
            return pending.get(index) != null && !pending.get(index).isEmpty();
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdrive.utils.Iterator#next()
         */
        public ExprProbeValue next() {
			ExprProbeValue epv = pending.get(index).removeFirst();
			return epv;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdrive.utils.Iterator#getDescriptor()
         */
        public String getDescriptor() {
            return "Expr Experiment Iterator:" + exptNames.get(index);
        } 
        
    }
}


