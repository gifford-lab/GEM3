package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;
import java.util.regex.*;
import java.io.*;

public class SOFTSample extends AbstractSOFTTableFile {
	
	private Map<String,Integer> idRefLines;
	private Integer valueIndex;
	private Boolean logValue; 

	public SOFTSample() { 
		super();
		idRefLines = new HashMap<String,Integer>();
		valueIndex = null;
		logValue = null;
	}
	
	public Set<String> allIDs() { 
		TreeSet<String> idSet = new TreeSet<String>();
		for(int i = 0; i < table.size(); i++) { 
			idSet.add(table.get(i)[0]);
		}
		return idSet;
	}
	
	public boolean hasID(String id) { 
		return idRefLines.containsKey(id);
	}
	
	public Integer findIDIndex(String id) { 
		return idRefLines.get(id); 
	}
	
	public Double getValue(int idx) {
		if(valueIndex == null) { 
			valueIndex = findValueHeaderIndex();
		}
		if(logValue == null) { 
			SOFTAttributes attrs = getTableAttributes();
			logValue = false;
			if(attrs.containsKey("VALUE")) { 
				String value = attrs.getCompleteValue("VALUE").toUpperCase();
				if(value.indexOf("LOG") != -1 
						|| value.indexOf("X SCORE") != -1
						|| value.indexOf("NORMALIZED") != -1
						|| value.indexOf("RMA") != -1) { 
					logValue = true;
				}
			}
		}
		try { 
			Double d = Double.parseDouble(table.get(idx)[valueIndex]);
			if(!logValue) { d = Math.log(d); }
			return d;
		} catch(NumberFormatException nfe) { 
			return null;
		} catch(ArrayIndexOutOfBoundsException aob) { 
			return null;
		}
	}
	
	public Double getValue(String id) {
		Integer idx = findIDIndex(id);
		if(idx != null) { 
			return getValue(findIDIndex(id));
		} else { 
			return null;
		}
	}
	
	public void printSampleLine(SOFTPlatform platform, PrintStream ps) {
		ps.print(String.format("%s", getAttributes().getCompleteValue("SAMPLE")));
		for(String id : platform.allIDs()) { 
			ps.print(String.format("\t%.3f", getValue(id)));
		}
		ps.println();
	}
	
	public int addToTable(String[] splitLine) { 
		int idx = super.addToTable(splitLine);
		String idRef = splitLine[0];
		if(idRefLines.containsKey(idRef)) { 
			String msg = String.format("Duplicate ID_REF: %s", idRef);
			throw new IllegalArgumentException(msg);
		}
		
		idRefLines.put(idRef, idx);
		return idx;
	}
}

	