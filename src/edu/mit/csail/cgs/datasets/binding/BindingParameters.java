package edu.mit.csail.cgs.datasets.binding;

import java.util.*;

public class BindingParameters {
	
	private Map<String,String> values;
	
	public BindingParameters() { 
		values = new HashMap<String,String>();
	}

	public BindingParameters(String paramString) { 
		values = new HashMap<String,String>();
		String[] array = paramString.split(";");
		for(int i = 0; i < array.length; i++) { 
			String[] pair = array[i].split("=");
			if(pair.length == 1) { 
				values.put(pair[0], "true");
			} else { 
				values.put(pair[0], pair[1]);
			}
		}
	}
	
	public boolean containsKey(String key) { return values.containsKey(key); }
	public String get(String key) { return values.get(key); }
	public void put(String key, String value) { values.put(key, value); }
	
	public double getDouble(String key, double defval) {
		try { 
			return values.containsKey(key) ? Double.parseDouble(values.get(key)) : defval;
		} catch(NumberFormatException nfe) { 
			return defval;
		}
	}
	
	public int getInt(String key, int defval) { 
		try { 
			return values.containsKey(key) ? Integer.parseInt(values.get(key)) : defval;
		} catch(NumberFormatException nfe) { 
			return defval;
		}
	}
	
	public String getString(String key, String defVal) { 
		return values.containsKey(key) ? values.get(key) : defVal;
	}
}
