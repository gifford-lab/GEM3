/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;

import edu.mit.csail.cgs.utils.models.Model;

public class BNValues {

	public Map<String,Object> values;
	
	public BNValues() { 
		values = new TreeMap<String,Object>();
	}
	
	public BNValues(Model m, BNVar... vars) {
		this();
		for(int i = 0; i < vars.length; i++) {
			if(values.containsKey(vars[i].getName())) { 
				throw new IllegalArgumentException(
						String.format("Duplicate name: %s", vars[i].getName()));
			}
			values.put(vars[i].getName(), vars[i].findValue(m));
		}
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		for(String k : values.keySet()) { 
			if(sb.length() > 1) { sb.append(" "); }
			sb.append(String.format("%s:%s", k, values.get(k).toString()));
		}
		sb.append("]");
		return sb.toString();
	}
	
	public BNValues(BNVar[] vars, Object[] vals) { 
		this();
		if(vars.length != vals.length) { throw new IllegalArgumentException(); }
		
		for(int i = 0; i < vars.length; i++) { 
			
			if(values.containsKey(vars[i].getName())) { 
				throw new IllegalArgumentException(
						String.format("Duplicate name: %s", vars[i].getName()));				
			}
			
			if(!vars[i].hasValue(vals[i])) {
				throw new IllegalArgumentException(vals[i].toString());
			}
			
			values.put(vars[i].getName(), vals[i]);
		}
	}
	
	public Object[] valueArray() {
		Object[] array = new Object[values.size()];
		int i = 0;
		for(String k : values.keySet()) { 
			array[i++] = values.get(k);
		}
		return array;
	}
	
	public String[] nameArray() { 
		String[] names = new String[values.size()];
		int i = 0; 
		for(String k : values.keySet()) { 
			names[i++] = k;
		}
		return names;
	}
	
	public int hashCode() { 
		int code = 17;
		for(String k : values.keySet()) { 
			code += k.hashCode(); code *= 37;
			code += values.get(k).hashCode(); code *= 37;
		}
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof BNValues)) { return false; }
		BNValues v = (BNValues)o;
		if(values.size() != v.values.size()) { return false; }
		for(String k : values.keySet()) { 
			if(!v.values.containsKey(k) || !v.values.get(k).equals(values.get(k))) { 
				return false; 
			}
		}
		return true;
	}
}
