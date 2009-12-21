package edu.mit.csail.cgs.utils.database;

import oracle.jdbc.OraclePreparedStatement;
import oracle.sql.*;
import java.sql.*;
import java.io.*;
import java.util.*;

public class ClobHandler {

	public static void setClob(Connection cxn, PreparedStatement ps, int index, String clobValue) throws SQLException { 
		if(DatabaseFactory.isOracle(cxn)) { 
			OraclePreparedStatement ops = (OraclePreparedStatement)ps;
			ops.setStringForClob(index, clobValue);
		} else { 
			ps.setString(index, clobValue);
		}
	}
}
