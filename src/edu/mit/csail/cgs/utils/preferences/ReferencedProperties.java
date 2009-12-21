package edu.mit.csail.cgs.utils.preferences;

import java.util.*;

public class ReferencedProperties {
	
	private ResourceBundle bundle;
	private Map<String,Object> props;

	public ReferencedProperties(String name) { 
		bundle = ResourceBundle.getBundle(name);
		props = new HashMap<String,Object>();
		Enumeration<String> keys = bundle.getKeys();
		while(keys.hasMoreElements()) { 
			String key = keys.nextElement();
			Object vals = dereferenceProperty(key);
			if(vals==null) { 
				System.err.println(String.format("Couldn't decode key: %s", key));
			} else { 
				props.put(key, vals);
			}
		}
	}
	
	public static String makeString(Object v) {
		if(v instanceof Vector) { 
			Vector vs = (Vector)v;
			StringBuilder sb = new StringBuilder();
			sb.append("[");
			for(int i = 0; i < vs.size(); i++) { 
				Object val = vs.get(i);
				if(i != 0) { sb.append(","); }
				sb.append(makeString(val));
			}
			sb.append("]");
			return sb.toString(); 
		} else { 
			return v.toString();
		}
	}

	public ResourceBundle getBundle() { return bundle; }
	
	public Object expand(String val) { 
		return dereferenceProperty(String.format("@%s", val));
	}
	
	private Object dereferenceProperty(String key) { 
		String[] ka = bundle.getString(key).split(",");
		Vector array = new Vector();
		for(int i = 0; i < ka.length; i++) { 
			array.add(dereferenceProperty(ka[i], new HashSet<String>()));
		}
		if(array.size() == 1) { 
			return array.get(0); 
		} else { 
			return array;
		}
	}
	
	private Object dereferenceProperty(String val, Set<String> visited) {
		if(!(val.startsWith("@"))) {
			return val;
		} else { 
			String key = val.substring(1, val.length());
			if(visited.contains(key)) { 
				return null;
			} else { 
				HashSet<String> newv = new HashSet<String>(visited);
				newv.add(key);
				String newval = bundle.getString(key);
				String[] nva = newval.split(",");
				Vector vec = new Vector();
				
				for(int i = 0; i < nva.length; i++) {
					if(nva[i].trim().length() > 0) { 
						Object obj = dereferenceProperty(nva[i], newv);
						vec.add(obj);
					}
				}
				
				if(nva.length == 1) {
					// we check nva.length instead of vec.size(), because we want to 
					// allow someone to specify a vector of length 1 with a syntax like
					// "foo=blah," in the properties file.  You'd get back a vector that
					// looks like "foo => [blah]", which would normally be "unwrapped".
					return vec.get(0);
				} else { 
					return vec;
				}
			}
		}
	}
	
	public Set<String> getKeys() { return props.keySet(); }
	
	public String getStringProperty(String key) { 
		return (String)getProperty(key);
	}
	
	public Vector getVectorProperty(String key) { 
		return (Vector)getProperty(key);
	}
	
	public Object getProperty(String key) { 
		return props.get(key); 
	}
	
	public boolean hasProperty(String key) { return props.containsKey(key); }
}
