/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.models.Model;

public class BNVar {

	private String name;
	private Map<Object,Integer> encoding;
	private Map<Integer,Object> decoding;
	
	public BNVar(String n, Set values) { 
		name = n;
		encoding = new HashMap<Object,Integer>();
		decoding = new TreeMap<Integer,Object>();
		
		int i = 0;
		for(Object v : values) { 
			encoding.put(v, i);
			decoding.put(i, v);
			i += 1;
		}
	}
	
	public Object findValue(Model m) { 
		try {
			Class mc = m.getClass();
			Field f = mc.getField(name);
			Class type = f.getType();
			if(Model.isSubclass(type, Integer.class) || Model.isSubclass(type, String.class)) {
				Object value = f.get(m);
				return value;
			}
			
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return null; 
	}
	
	public void setValue(Model m, Object value) { 
		try {
			Class mc = m.getClass();
			Field f = mc.getField(name);
			Class type = f.getType();
			if(Model.isSubclass(type, Integer.class) || Model.isSubclass(type, String.class)) {
				f.set(m, value);
			}
			
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
	}
	
	public String getName() { return name; }
	public int size() { return encoding.size(); }
	
	public Integer encode(Object value) { 
		return encoding.get(value);
	}
	
	public Object decode(Integer i) { 
		return decoding.get(i);
	}
	
	public boolean hasValue(Object v) { return encoding.containsKey(v); }
	
	public String toString() { return name; }
	
	public int hashCode() { 
		return name.hashCode(); 
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof BNVar)) { return false; }
		BNVar b = (BNVar)o;
		return b.name.equals(name);
	}
}
