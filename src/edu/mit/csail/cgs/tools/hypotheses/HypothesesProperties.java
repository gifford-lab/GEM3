package edu.mit.csail.cgs.tools.hypotheses;

import java.util.*;

public class HypothesesProperties {
	
	private static ResourceBundle bundle;
	private static HashSet<String> bundleKeys;
	
	static { 
		bundle = ResourceBundle.getBundle("edu.mit.csail.cgs.tools.hypotheses.hypotheses");
		bundleKeys = new HashSet<String>();
		Enumeration<String> keys = bundle.getKeys();
		while(keys.hasMoreElements()) { 
			bundleKeys.add(keys.nextElement());
		}
	}
	
	public static Set<String> getKeys() { 
		return bundleKeys;
	}
	
	public static boolean hasKey(String k) { 
		return bundleKeys.contains(k); 
	}
	
	public static String getValue(String k) { 
		return bundle.getString(k);
	}
}
