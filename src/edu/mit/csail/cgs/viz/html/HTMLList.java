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
public class HTMLList implements HTMLElmt {
	
	private String listType;
	private LinkedList<HTMLElmt> elmtList;
	
	public HTMLList(boolean ordered) { 
		if(ordered) { 
			listType = "ol";
		} else { 
			listType = "ul";
		}
		
		elmtList = new LinkedList<HTMLElmt>();
	}
	
	public HTMLList(boolean ordered, Collection<HTMLElmt> elmts) { 
		if(ordered) { 
			listType = "ol";
		} else { 
			listType = "ul";
		}
		
		elmtList = new LinkedList<HTMLElmt>(elmts);		
	}
	
	public void addElmt(HTMLElmt elmt) { elmtList.addLast(elmt); }
	public void clear() { elmtList.clear(); }
	public int size() { return elmtList.size(); }
	public String getListType() { return listType; }

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.html.HTMLElmt#print(java.io.PrintStream)
	 */
	public void print(PrintStream ps) {
		ps.println("<" + listType + ">");
		for(HTMLElmt elmt : elmtList) { 
			ps.print("<li>");
			elmt.print(ps);
			ps.print("</li>");
		}
		ps.println("</" + listType + ">");
	}

}
