package edu.mit.csail.cgs.tools.sql;

public interface SQLType {
	
	public static final String ORACLE_TYPE = "oracle";
	public static final String MYSQL_TYPE = "mysql";
	
	public String translateToDB(String db);
	
	public static class StringType implements SQLType { 
		private int size;
		
		public StringType(int s) { 
			size = s;
		}

		public String translateToDB(String db) { 
			if(db.equals(ORACLE_TYPE)) { 
				return "varchar2(" + size + ")";
			}
			
			if(db.equals(MYSQL_TYPE)) { 
				return "varchar(" + size + ")";
			}
			
			throw new UnsupportedDatabaseException(db);
		}
	}
	
	public static class IntegerType implements SQLType {
		
		private int size;
		
		public IntegerType(int s) { size = s; }
		
		public String translateToDB(String db) { 
			if(db.equals(ORACLE_TYPE)) { 
				return "number(" + size + ")";
			}
			
			if(db.equals(MYSQL_TYPE)) { 
				return "int(" + size + ")";
			}
			
			throw new UnsupportedDatabaseException(db);
		}
	}
}
