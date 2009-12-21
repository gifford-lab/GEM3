/**
 * 
 */
package edu.mit.csail.cgs.viz.html;

import java.io.PrintStream;
import java.util.Collection;
import java.util.LinkedList;

/**
 * @author Timothy Danford
 *
 */
public class HTMLTable implements HTMLElmt {
	
	private LinkedList<TR> rowList;
	
	public HTMLTable() { rowList = new LinkedList<TR>(); }
	
	public HTMLTable(LinkedList<TR> rows) { rowList = new LinkedList<TR>(rows); }
	
	public HTMLTable(HTMLElmt[][] array) { 
		rowList =new LinkedList<TR>();
		for(int i = 0; i < array.length; i++) { 
			rowList.addLast(new TR(array[i]));
		}
	}
	
	public LinkedList<LinkedList<HTMLElmt>> getElmts() { 
		LinkedList<LinkedList<HTMLElmt>> elmts = new LinkedList<LinkedList<HTMLElmt>>();
		for(TR tr : rowList) { 
			elmts.addLast(tr.getElmtList());
		}
		return elmts;
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.cgs.viz.html.HTMLElmt#print(java.io.PrintStream)
	 */
	public void print(PrintStream ps) {
		ps.print("<table>");
		for(TR row : rowList) { 
			row.print(ps);
		}
		ps.println("</table>");
	}

}

class TR { 
	private LinkedList<TD> colList;
	
	public TR() { 
		colList = new LinkedList<TD>(); 
	}
	
	public TR(Collection<HTMLElmt> elmtList) { 
		colList = new LinkedList<TD>();
		for(HTMLElmt elmt : elmtList) { 
			colList.addLast(new TD(elmt)); 
		}
	}
	
	public TR(HTMLElmt[] array) { 
		colList = new LinkedList<TD>();
		for(int i = 0; i < array.length; i++) { 
			colList.addLast(new TD(array[i]));
		}
	}
	
	public void addElmt(HTMLElmt elmt) { colList.addLast(new TD(elmt)); }
	
	public LinkedList<HTMLElmt> getElmtList() { 
		LinkedList<HTMLElmt> elmts =new LinkedList<HTMLElmt>();
		for(TD td : colList) { elmts.addLast(td.getElmt()); }
		return elmts;
	}
	
	public void print(PrintStream ps) { 
		ps.print("<tr>");
		for(TD td : colList) { 
			td.print(ps);
		}
		ps.println("</tr>");
	}
}

class TD { 
	private HTMLElmt elmt;
	
	public TD(HTMLElmt elmt) { this.elmt = elmt; }
	public HTMLElmt getElmt() { return elmt; }
	public void setElmt(HTMLElmt elmt) { this.elmt = elmt; }
	
	public void print(PrintStream ps) { 
		ps.print("<td>");
		elmt.print(ps);
		ps.print("</td>");
	}
}