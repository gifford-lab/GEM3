/*
 * Created on Jan 7, 2008
 */
package edu.mit.csail.cgs.utils.parsing.alignment;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.Predicate;

public class ElandFile implements Iterator<ElandHit> {
    
    public static void main(String[] args) { 
        File f = new File(args[0]);
        try {
            ElandFile eland = new ElandFile(f);
            int c = 0;
            while(eland.hasNext()) { 
                ElandHit hit = eland.next();
                c += 1;
            }
            
            System.out.println(String.format("# Eland Hits: %d", c));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private Set<ElandHit.Code> codes;
    private ElandHit nextHit;
    private BufferedReader br;
    private Predicate<ElandHit> predicate;
    private int lines;

    public ElandFile(File f, Collection<ElandHit.Code> cds, Predicate<ElandHit> pred) throws IOException { 
        codes = new HashSet<ElandHit.Code>(cds);
        nextHit = null;
        lines = 0;
        br = new BufferedReader(new FileReader(f));
        predicate = pred;
        findNextHit();
    }


    public ElandFile(File f, Collection<ElandHit.Code> cds) throws IOException { 
        codes = new HashSet<ElandHit.Code>(cds);
        nextHit = null;
        lines = 0;
        br = new BufferedReader(new FileReader(f));
        predicate = null;
        findNextHit();
    }

    public ElandFile(File f) throws IOException { 
        codes = new HashSet<ElandHit.Code>();
        codes.add(ElandHit.Code.U0);
        codes.add(ElandHit.Code.U1);
        //codes.add(ElandHit.Code.U2);
        
        nextHit = null;
        lines = 0;
        predicate = null;
        br = new BufferedReader(new FileReader(f));
        findNextHit();
    }
    
    public boolean acceptsHit(ElandHit hit) { 
    	return codes.contains(hit.getCode()) && 
    		(predicate == null || predicate.accepts(hit));
    }

    private void findNextHit() {
        
        String line;
        boolean error = false;
        do {
        	do { 
	            nextHit = null;
	            error = false;
	            try {
	                line = br.readLine();
	                if(line != null) { 
	                    try {
	                        nextHit = new ElandHit(lines, line);
	                        if(nextHit==null || nextHit.getName() == null || nextHit.getQuerySequence() == null) { 
	                            error = true;
	                        }
	                    } catch (ElandParsingException e) {
	                        error = true;
	                    }
	                    if(error) { 
	                        System.err.println(String.format(
	                                "Line #%d had an error; entry was ignored.", lines));
	                    }
	                    lines++;
	                }
	            } catch (IOException e) {
	            	error=true;
	                //e.printStackTrace();
	            }
        	}while(error);
        } while(nextHit != null && (error || !acceptsHit(nextHit)));
        
        if(nextHit==null) { 
            try {
                br.close();
            } catch (IOException e) {
                e.printStackTrace();
            } 
        }
    }

    public boolean hasNext() {
        return nextHit != null;
    }

    public ElandHit next() {
        ElandHit nexthit = nextHit;
        findNextHit();
        return nexthit;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

}
