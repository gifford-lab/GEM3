package edu.mit.csail.cgs.tools.sql;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.regex.*;

public class SQLParser {
	
	public static void main(String[] args) {
		String fname = "../gse-db-schemas/annotations.oracle";
		File f = new File(fname);
		try {
			SQLParser parser = new SQLParser(f);
			parser.printStatements();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private Vector<String> statements;

	public SQLParser(File f) throws IOException { 
		BufferedReader br = new BufferedReader(new FileReader(f));
		StringBuilder sb = new StringBuilder();
		String line = null;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			if(line.length() > 0 && !line.startsWith("--")) { 
				sb.append(line);
			}
		}
		br.close();
		
		String[] stmts = sb.toString().split(";");
		statements = new Vector<String>();
		for(int i = 0; i < stmts.length; i++) { 
			statements.add(stmts[i]);
		}
	}
	
	public static Object parseStatement(String stmt) { 
		stmt = stmt.toLowerCase();
		Pattern p = Pattern.compile("^create\\s+table\\s+.*");
		Matcher m = p.matcher(stmt);
		if(m.matches()) { return new SQLCreateTableStatement(stmt); }
		
		return null;
	}
	
	public Vector<SQLCreateTableStatement> getTableCreationStatements() { 
		Vector<SQLCreateTableStatement> stmts = new Vector<SQLCreateTableStatement>();
		for(String stmt : statements) { 
			Object parsed = parseStatement(stmt);
			if(parsed != null && (parsed instanceof SQLCreateTableStatement)) { 
				stmts.add((SQLCreateTableStatement)parsed);
			}
		}
		return stmts;
	}
	
	public void printStatements() {
		int i = 0;
		for(String stmt : statements) {
			Object parsed = parseStatement(stmt);
			if(parsed == null) { 
				System.out.println(i + ": \"" + stmt +"\"");
			} else { 
				System.out.println(i + ":\n" + parsed.toString());
			}
			i++;
		}
	}
}
