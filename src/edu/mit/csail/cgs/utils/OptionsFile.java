package edu.mit.csail.cgs.utils;

import java.io.*;
import java.util.*;

public class OptionsFile {
	
	private Map<String,Vector<String>> options;

	public OptionsFile(File f) throws IOException { 
		options = new HashMap<String,Vector<String>>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		while((line = br.readLine()) != null) { 
		    if(line.indexOf("=") != -1) { 
		        String[] array = line.split("=");
		        String key = array[0].trim();
		        String values = array[1];
		        Vector<String> vals = new Vector<String>();
		        String[] varray = values.split(",");
		        for(int i = 0; i < varray.length; i++) { vals.add(varray[i].trim()); }
		        options.put(key, vals);
		    }
		}
		
		br.close();
	}

	public int size() { return options.size(); }
	public String getValue(String k) { return options.get(k).get(0); }
	public String getValue(String k, int i) { return options.get(k).get(i); }
	public Vector<String> getAllValues(String k) { return options.get(k); }
	public boolean hasValue(String k) { return options.containsKey(k); }
	public Set<String> getKeys() { return options.keySet(); }
}
