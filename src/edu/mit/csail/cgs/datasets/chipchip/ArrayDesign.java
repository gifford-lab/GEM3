/**
 * 
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author tdanford
 */
public class ArrayDesign {
    
    public static void main(String[] args) { 
        try {
            Collection<ArrayDesign> values = loadDBArrayDesigns();
            for(ArrayDesign v : values) { System.out.println(v.toString()); }
            
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }        
    }
	
	public static Collection<ArrayDesign> loadDBArrayDesigns()
		throws SQLException, UnknownRoleException { 
		
		LinkedList<ArrayDesign> cells = new LinkedList<ArrayDesign>();
		java.sql.Connection c = 
			DatabaseFactory.getConnection(ExptLocator.dbRole);
		Statement s = c.createStatement();
		ResultSet rs = s.executeQuery("select id, name, genome from arraydesign");
		
		while(rs.next()) { 
			ArrayDesign cell = new ArrayDesign(rs);
			cells.addLast(cell);
		}
		
		rs.close();
		s.close();
		DatabaseFactory.freeConnection(c);
		return cells;
	}

	private int dbid;
	private String name;
    private Genome genome;
    
    public ArrayDesign(ResultSet rs) throws SQLException { 
        dbid = rs.getInt(1);
        name = rs.getString(2);
        int genomeID = rs.getInt(3);
        try {
            genome = Organism.findGenome(genomeID);
        } catch (NotFoundException e) {
            throw new IllegalArgumentException("Unknown Genome: " + genomeID + " was in ArrayDesign table!", e);
        }
    }
	
	public String getName() { return name; }
	public int getDBID() { return dbid; }
    public Genome getGenome() { return genome; }
    
    public String toString() { return name + " (#" + dbid + ") --> " + genome.getVersion(); }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ArrayDesign)) { return false; }
		ArrayDesign c = (ArrayDesign)o;
		if(dbid != c.dbid) { return false; }
		return true;
	}
}
