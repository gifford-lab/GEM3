/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.io.PrintStream;
import java.util.*;

/**
 * @author Timothy Danford
 *
 */
public class HTMLSequence implements HTMLElmt {
	
	private LinkedList<HTMLElmt> elmts;
	
	public HTMLSequence(Collection<HTMLElmt> elmts) { 
		this.elmts =new LinkedList<HTMLElmt>(elmts);
	}
	
	public HTMLSequence() { elmts =new LinkedList<HTMLElmt>(); }
	
	public void addElmt(HTMLElmt elmt) { elmts.addLast(elmt); }
	public void clear() { elmts.clear(); }
	
	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.html.HTMLElmt#print(java.io.PrintStream)
	 */
	public void print(PrintStream ps) {
		for(HTMLElmt elmt : elmts) { elmt.print(ps); }
	}

}
