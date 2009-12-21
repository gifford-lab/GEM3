/*
 * Author: tdanford
 * Date: Nov 6, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;

/**
 * The contents of this class are just helper functions -- to check if all required fields
 * are present, and to print errors if they are not.  
 * 
 * In general, all a subclass needs to do is to extend ArgumentModel, provide a reasonable
 * default constructor if necessary, and then make sure the Class object of the subclass 
 * is passed into the appropriate Arguments object.  Arguments.java does all the parsing for us.  
 * 
 * @author tdanford
 */
public class ArgumentModel extends Model {
	
	private ModelFieldAnalysis analysis;
	
	public ArgumentModel() { 
		analysis = new ModelFieldAnalysis(getClass());
	}
	
	public List<String> findRequiredFields() {
		ArrayList<String> reqlist = new ArrayList<String>();

		Field f = analysis.findStaticField("required");
		if(f != null && Model.isSubclass(f.getType(), String[].class)) { 
			try {
				String[] reqarray = (String[]) f.get(this);
				for(int i = 0; reqarray != null && i < reqarray.length; i++) { 
					reqlist.add(reqarray[i]);
				}
				
			} catch (IllegalAccessException e) {
				// do nothing.
			}
		}
		
		return reqlist;
	}
	
	public boolean checkArgs() {
		List<String> args = findRequiredFields();
		
		for(String arg : args) { 
			Field f = analysis.findField(arg);
			if(f == null) { return false; }
			try {
				Object value = f.get(this);
				if(value == null) { return false; }
				
			} catch (IllegalAccessException e) {
				return false;
			}
		}
		return true;
	}
	
	public String getArgErrors() { 
		StringBuilder sb = new StringBuilder();
		List<String> args = findRequiredFields();
		sb.append("Missing arguments:");
		
		argloop: for(String arg : args) { 
			Field f = analysis.findField(arg);
			if(f == null) { 
				sb.append(" " + arg);
				continue argloop;
			}
			
			try {
				Object value = f.get(this);
				if(value == null) { 
					sb.append(" " + arg);
					continue argloop;
				}
				
			} catch (IllegalAccessException e) {
				sb.append(" " + arg); 
				continue argloop;
			}
		}

		return sb.toString();
	}
}
