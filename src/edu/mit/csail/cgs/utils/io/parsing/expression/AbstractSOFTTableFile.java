package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;
import java.util.regex.*;
import java.io.*;

public class AbstractSOFTTableFile extends AbstractSOFTFile {
	
	public static String splitToken = "\t";

	protected SOFTAttributes tableAttrs;
	protected Vector<String> header;
	protected Vector<String[]> table;
	
	public AbstractSOFTTableFile() { 
		super();
		tableAttrs = new SOFTAttributes();
		header = null;
		table = new Vector<String[]>();
	}
	
	public int findHeaderIndex(String h) { return header.indexOf(h); }
	public int findValueHeaderIndex() { return findHeaderIndex("VALUE"); }
	
	public void printHeaderLine(PrintStream ps) {
		for(int i = 0; i < header.size(); i++) { 
			ps.print(String.format("%s%s", i > 0 ? "\t" : "", header.get(i)));
		}
	}
	
	public void printTableLine(PrintStream ps, int idx) { 
		String[] line = table.get(idx);
		for(int i = 0; i < line.length; i++) { 
			ps.print(String.format("%s%s", i > 0 ? "\t" : "", line[i]));
		}
	}
	
	public SOFTAttributes getTableAttributes() { return tableAttrs; }
	
	public void setHeader(String line) {
		if(header != null) { 
			String msg = String.format("Tried to set header values twice (%s)",
					header.toString());
			throw new IllegalArgumentException(msg);
		}
		
		String[] a = line.split(splitToken);
		header = new Vector<String>();
		for(int i = 0; i < a.length; i++) {
			header.add(a[i]);
		}
	}
	
	public int addToTable(String line) {
		String[] a = line.split(splitToken);
		try { 
			return addToTable(a);
		} catch(IllegalArgumentException e) { 
			String msg = String.format("\"%s\" -> '%s'", 
					e.getMessage(), line);
			throw new IllegalArgumentException(msg, e);
		}
	}
	
	public int addToTable(String[] splitLine) { 
		if(header==null) {  
			throw new IllegalArgumentException();
		}
		if(splitLine.length > header.size()) {
			for(int i = 0; i < header.size(); i++) { 
				System.err.print(i + ":" + header.get(i) + " ");
			}
			System.err.println();
			for(int i = 0; i < splitLine.length; i++) { 
				System.err.print(i + ":" + splitLine[i] + " ");
			}
			System.err.println();
			throw new IllegalArgumentException(String.format("%d != %d", 
					splitLine.length, header.size()));
		}

		table.add(splitLine);
		return table.size()-1; 
	}
	
	public String[] getTableRow(int i) { return table.get(i); }
	public int getTableSize() { return table.size(); }
	public Vector<String> getHeader() { return header; }
}
