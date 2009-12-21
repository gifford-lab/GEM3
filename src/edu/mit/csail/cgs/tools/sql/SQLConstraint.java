package edu.mit.csail.cgs.tools.sql;

public interface SQLConstraint {

	public String translate(String db);
	
	public static class ReferencesConstraint implements SQLConstraint {
		
		private String tableName;
		private String fieldName;
		
		public ReferencesConstraint(String tn, String fn) { 
			tableName = tn;
			fieldName = fn;
		}
		
		public String translate(String db) { 
			return "references " + tableName + "(" + fieldName + ")";
		}
	}
}
