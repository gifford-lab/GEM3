package edu.mit.csail.cgs.tools.sql;

import java.util.*;
import java.util.regex.*;

public class SQLCreateTableStatement {
	
	private String name;
	private Vector<TableField> fields;
	private String modifiers;
	
	public SQLCreateTableStatement(String stmt) { 
		stmt = stmt.toLowerCase();
		Pattern p = Pattern.compile("create table ([\\_\\w+\\d+]+)\\s*\\((.*)\\)\\s*(.*)");
		Matcher m = p.matcher(stmt);
		if(m.matches()) { 
			name = m.group(1);
			modifiers = m.group(3);
			String fieldString = m.group(2);
			fields = new Vector<TableField>();
			
			int comma = -1, sparen = -1, eparen = -1;
			comma = fieldString.indexOf(",");
			sparen = fieldString.indexOf("(");
			eparen = fieldString.indexOf(")");
			
			while(comma != -1) { 
				//System.out.println(comma + "," + sparen + "," + eparen + " \"" + fieldString + "\"");

				if(sparen != -1 && sparen < comma && eparen > comma) {  
					fieldString = fieldString.substring(eparen+1, fieldString.length());
				} else {
					if(comma > 0) { 
						String fieldstr = fieldString.substring(0, comma);
						if(!fieldstr.startsWith("constraint")) { 
							fields.add(new TableField(fieldstr));
						} else { 
							modifiers += ", " + fieldstr;
						}
					}
					fieldString = fieldString.substring(comma+1, fieldString.length());
				}

				comma = fieldString.indexOf(",");
				sparen = fieldString.indexOf("(");
				eparen = fieldString.indexOf(")");
			}
			
			if(fieldString.length() > 0) { 
				fields.add(new TableField(fieldString));
			}
			
		} else { 
			throw new IllegalArgumentException("Illegal Table Statement: \"" + stmt + "\"");
		}
	}
	
	public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("Table \"" + name + "\" --> (" + modifiers + ")\n");
		for(TableField tf : fields) { sb.append("\t" + tf + "\n"); }
		return sb.toString();
	}
}

class TableField { 
	private String name, type, modifiers;
	
	public TableField(String field) { 
		Pattern p = Pattern.compile("(\\w+)\\s+([\\w]+(?:\\([\\s\\d]+\\))?)(.*)");
		Matcher m = p.matcher(field);
		if(m.matches()) { 
			name = m.group(1);
			type = m.group(2);
			modifiers = m.group(3);
			
			if(type.equals("constraint")) { 
				modifiers = type + " " + modifiers;
				type = "";
			}
		} else {
			throw new IllegalArgumentException("Unknown Field \"" + field + "\"");
		}
	}
	
	public TableField(String n, String t, String m) { 
		name = n;
		type = t;
		modifiers = m;
	}
	
	public String toString() { 
		return "\"" + name + "\" (" + type + ") --> " + modifiers;
	}
}
