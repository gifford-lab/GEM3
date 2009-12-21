/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.util.*;
import java.io.*;

/**
 * @author Timothy Danford
 */
public class HTMLPage {
	
	private String title;
	private LinkedList<HTMLElmt> elmtList;
	
	public HTMLPage(String title) { 
		this.title = title;
	}
	
	public void addElmt(HTMLElmt elmt) { 
		elmtList.addLast(elmt);
	}
	
	public void printPage(PrintStream ps) {
		ps.println("<html>");
		ps.println("<head><title>" + title + "</title></head>");
		ps.println("<body>");
		for(HTMLElmt elmt : elmtList) { 
			elmt.print(ps);
		}
		ps.println("</body>\n</html>");
	}
}
