/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.io.PrintStream;

/**
 * @author Timothy Danford
 *
 */
public class HTMLBreak implements HTMLElmt {
	
	public HTMLBreak() {}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.html.HTMLElmt#print(java.io.PrintStream)
	 */
	public void print(PrintStream ps) {
		ps.print("<br>");
	}

}
