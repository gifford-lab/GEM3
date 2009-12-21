/*
 * Created on May 23, 2005
 */
package edu.mit.csail.cgs.utils.parsing;

import java.util.*;

import java.sql.*;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public class ParseDBExpr implements ParseExpr {
    
    private java.sql.Connection cxn;
    private Statement stmt;
    
    public ParseDBExpr() { 
        try {
            java.sql.Connection cxn = DatabaseFactory.getConnection("psrg");
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role psrg",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to psrg",ex);
        }
        try { 
            stmt = cxn.createStatement();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }
    
    public void close() {
        try { 
            stmt.close();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
        DatabaseFactory.freeConnection(cxn);
    }

    public String getExptName(int i) {
        try { 
            ResultSet rs = stmt.executeQuery("select name from expr_experiment where id=" + i);
            if(rs.next()) { 
                String n = rs.getString(1);
                rs.close();
                return n;
            } else { 
                rs.close();
                throw new IllegalArgumentException("No expr_experiment: " + i);
            }
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

    public Map<String, Double> getWholeExpt(int i) {
        Map<String,Double> vals = new HashMap<String,Double>();
        try { 
            ResultSet rs = stmt.executeQuery("select probe_name, value from expr where expt_index = " + i +
				" and value is not null");
            while(rs.next()) { 
                vals.put(rs.getString(1), rs.getDouble(2));
            }
            rs.close();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
        
        return vals;
    }

    public boolean exptHasORF(int i, String orf) {
        try { 
            ResultSet rs = stmt.executeQuery("select count(*) from expr where probe_name = '" + orf + "' and expt_index=" + i + " and value is not null");
            int count = 0;
            if(rs.next()) { 
                count = rs.getInt(1);
            }
            rs.close();
            return count > 0;
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

	public Double[] getValues(int start, int end, String orf) { 
		Double[] array = new Double[end-start+1];
		try { 
			for(int i = 0; i < array.length; i++) { array[i] = null; }
			ResultSet rs = stmt.executeQuery("select ee.id, e.value from expr e, expr_experiment ee " + 
				"where e.expt_index=ee.id and e.probe_name='" + orf + "' and ee.id >= " + start + 
				" and ee.id <= " + end + " and e.value is not null");
			while(rs.next()) { 
				int id = rs.getInt(1);
				double value = rs.getDouble(2);
				array[id-start] = new Double(value);
			}
			rs.close();
			return array;
		} catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
		}
	}

    public double getValue(int i, String orf) {
        try { 
            ResultSet rs = stmt.executeQuery("select value from expr where probe_name = '" + orf + "' and expt_index=" + i + " and value is not null");
            if(rs.next()) { 
                double v= rs.getDouble(1);
                rs.close();
                return v;
            } else { 
                rs.close();
                throw new IllegalArgumentException("No such ORF " + orf + " in expt " + i);
            }
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

    public int numExpts() {
        try { 
            ResultSet rs = stmt.executeQuery("select count(*) from expr_experiment");
            int count = 0;
            if(rs.next()) { 
                count = rs.getInt(1);
            }
            rs.close();
            return count;
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

    public int numORFs(int i) {
        try { 
            ResultSet rs = stmt.executeQuery("select count(*) from expr where expt_index=" + i + " and value is not null");
            int count = 0;
            if(rs.next()) { 
                count = rs.getInt(1);
            }
            rs.close();
            return count;
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

    public double getMax(int i) {
        try { 
            ResultSet rs = stmt.executeQuery("select max(value) from expr where expt_index=" + i + " and value is not null");
            double extremum = 0;
            if(rs.next()) { 
                extremum = rs.getDouble(1);
            }
            rs.close();
            return extremum;
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
    }

    public double getMin(int i) {
        try { 
            ResultSet rs = stmt.executeQuery("select min(value) from expr where expt_index=" + i + " and value is not null");
            double extremum = 0;
            if(rs.next()) { 
                extremum = rs.getDouble(1);
            }
            rs.close();
            return extremum;
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new RuntimeException(se.getMessage());
        }
   }
}
