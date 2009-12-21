/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.io.PrintStream;

/**
 * @author Timothy Danford
 *
 */
public class HTMLText implements HTMLElmt {
	
	private StringBuffer buffer;
	
	public HTMLText(String s) { 
		buffer = new StringBuffer(s);
	}
	
	public void addText(String s) { 
		buffer.append(s);
	}
	
	public void clear() { buffer = new StringBuffer(); }

	public void print(PrintStream ps) {
		ps.println(buffer.toString());
	}

}
