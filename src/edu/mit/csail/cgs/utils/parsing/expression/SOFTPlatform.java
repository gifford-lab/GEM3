package edu.mit.csail.cgs.utils.parsing.expression;

import java.util.*;
import java.util.regex.*;
import java.io.*;

public class SOFTPlatform extends AbstractSOFTTableFile {
	
	private Map<String,Integer> idMap;

	public SOFTPlatform() { 
		super();
		idMap = new HashMap<String,Integer>();
	}
	
	public boolean hasID(String id) { 
		return idMap.containsKey(id);
	}
	
	public Set<String> allIDs() { 
		TreeSet<String> idSet = new TreeSet<String>();
		for(int i = 0; i < table.size(); i++) { 
			idSet.add(table.get(i)[0]);
		}
		return idSet;
	}
	
	public void printSampleHeaderLine(PrintStream ps) {
		ps.print("Sample");
		for(String id : allIDs()) { 
			ps.print(String.format("\t%s", id));
		}
		ps.println();
	}

	public int findIDIndex(String id) { return idMap.get(id); }
	
	public int addToTable(String[] splitLine) { 
		int idx = super.addToTable(splitLine);
		String id = splitLine[0];
		if(idMap.containsKey(id)) { 
			String msg = String.format("Duplicate ID: %s", id);
			throw new IllegalArgumentException(msg);
		}
		
		idMap.put(id, idx);
		return idx;
	}
	
}
