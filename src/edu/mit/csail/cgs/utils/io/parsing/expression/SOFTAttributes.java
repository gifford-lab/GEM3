package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;
import java.util.regex.*;

public class SOFTAttributes { 
	private static Pattern linePattern = 
		Pattern.compile("!([^=]+)\\s*=\\s*(.*)");
	
	private Map<String,Vector<String>> values;
	
	public SOFTAttributes() { 
		values = new LinkedHashMap<String,Vector<String>>();
	}
	
	public boolean containsKey(String k) { 
		return values.containsKey(k);
	}
	
	public Collection<String> getAllValues(String k) { 
		return values.get(k);
	}
	
	public String getFirstValue(String k) { 
		return values.get(k).get(0);
	}
	
	public String getCompleteValue(String k) { 
		StringBuilder sb = new StringBuilder();
		Vector<String> lines = values.get(k);
		
		for(int i = 0; i < lines.size(); i++) { 
			sb.append(i > 0 ? " " : "");
			sb.append(lines.get(i));
		}
		
		return sb.toString();
	}
	
	public void addKeyValue(String key, String value) { 
		if(!values.containsKey(key)) { 
			values.put(key, new Vector<String>());
		}
		
		values.get(key).add(value);		
	}
	
	public void addLine(String line) { 
		Matcher m = linePattern.matcher(line);
		if(!m.matches()) { 
			throw new IllegalArgumentException(line);
		}
		
		String key = m.group(1);
		String value = m.group(2);
		
		addKeyValue(key, value);
	}
}
