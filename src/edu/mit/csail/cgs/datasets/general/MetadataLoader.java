package edu.mit.csail.cgs.datasets.general;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

/* to make MetadataLoader thread-safe, many of its method have an internal block
   that is synchronized on some PreparedStatement.  To prevent deadlocks, we want a 
   directed acyclic graph of which Loader uses which.

   Please add entries below for other loaders that also synchronize on themselves or fields
   to ensure that we don't create cycles (indentation below an entry indicates "uses").

   MetadataLoader
       ChipChipMetadataLoader
       ExpressionMetadataLoader
       ExpressionLoader
           ProbeMappingLoader
   TimeSeriesLoader
       ExpressionLoader
           ProbeMappingLoader
   ChipPetLoader
   BindingScanLoader
   LocatorLoader
   OrthologyLoader (uses itself through OrthologyPair constructor)
   TextFunctionLoader
   function.DatabaseFunctionLoader
   WeightMatrixLoader
   DistributionLoader
   AlignmentLoader
       ConservationLoader (uses an AlignmentLoader)
   

*/

public class MetadataLoader implements edu.mit.csail.cgs.utils.Closeable {
    
    public static final String role = "core";

    private Map<String,Cells> cellsNames;
    private Map<String,Condition> condNames;
    private Map<String,Factor> factorNames;
	
    private Map<Integer,Cells> cellsIDs;
    private Map<Integer,Condition> condIDs;
    private Map<Integer,Factor> factorIDs;
	
    private java.sql.Connection cxn;
    
    private PreparedStatement loadCells, loadCond, loadFactor;
    private PreparedStatement loadCellsByName, loadCondByName, loadFactorByName;
	
    public MetadataLoader() throws SQLException { 
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
        
        loadCells = cxn.prepareStatement("select id, name from cells where id=?");
        loadCond = cxn.prepareStatement("select id, name from conditions where id=?");
        loadFactor = cxn.prepareStatement("select id, name from factors where id=?");

        loadCellsByName = cxn.prepareStatement("select id, name from cells where name=?");
        loadCondByName = cxn.prepareStatement("select id, name from conditions where name=?");
        loadFactorByName = cxn.prepareStatement("select id, name from factors where name=?");
        
        cellsNames = new HashMap<String,Cells>();
        cellsIDs = new HashMap<Integer,Cells>();
        condNames = new HashMap<String,Condition>();
        condIDs = new HashMap<Integer,Condition>();
        factorNames = new HashMap<String,Factor>();
        factorIDs = new HashMap<Integer,Factor>();
    }
	
    public boolean isClosed() { return cxn==null; }

    public void close() { 
        try {
            loadCells.close(); loadCells = null;
            loadCond.close();  loadCond = null;
            loadFactor.close(); loadFactor = null;
            
            loadCellsByName.close();  loadCellsByName = null;
            loadCondByName.close(); loadCondByName = null;
            loadFactorByName.close(); loadFactorByName = null;
            
        } catch (SQLException e) {
            e.printStackTrace();
        }
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }
    
    public java.sql.Connection getConnection() { return cxn; }
    
    public Cells getCells(String name) throws SQLException { 
        synchronized (loadCellsByName) {
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                Cells c = new Cells(rs);
                rs.close();
                
                if(!cellsIDs.containsKey(c.getDBID())) { 
                    cellsIDs.put(c.getDBID(), c);
                    cellsNames.put(c.getName(), c);
                }
                
                return c;
            }
            rs.close();
        }
        int id = insertCells(name);
        return loadCells(id);
    }
    
    public Cells findCells(String name) throws SQLException { 
        synchronized (loadCellsByName) {
            loadCellsByName.setString(1, name);
            ResultSet rs = loadCellsByName.executeQuery();
            
            if(rs.next()) { 
                Cells c = new Cells(rs);
                rs.close();
                
                if(!cellsIDs.containsKey(c.getDBID())) { 
                    cellsIDs.put(c.getDBID(), c);
                    cellsNames.put(c.getName(), c);
            }
                
                return c;
            }            
            rs.close();
            return null;
        }
    }
    
    public Condition getCondition(String name) throws SQLException { 
        synchronized (loadCondByName) {
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                Condition c = new Condition(rs);
                rs.close();
                
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
                
                return c;
            }
            rs.close();
        }
        int id = insertCondition(name);
        return loadCondition(id);
    }        

    public Condition findCondition(String name) throws SQLException { 
        synchronized (loadCondByName) {
            loadCondByName.setString(1, name);
            ResultSet rs = loadCondByName.executeQuery();
            
            if(rs.next()) { 
                Condition c = new Condition(rs);
                rs.close();
                
                if(!condIDs.containsKey(c.getDBID())) { 
                    condIDs.put(c.getDBID(), c);
                    condNames.put(c.getName(), c);
                }
                
                return c;
            }            
            rs.close();
            return null;
        }
    }
    
    public Factor findFactor(String name) throws SQLException { 
        synchronized (loadFactorByName) {
            loadFactorByName.setString(1, name);
            ResultSet rs = loadFactorByName.executeQuery();
            
            if(rs.next()) { 
                Factor c = new Factor(rs);
                rs.close();
                
                if(!factorIDs.containsKey(c.getDBID())) { 
                    factorIDs.put(c.getDBID(), c);
                    factorNames.put(c.getName(), c);
                }
                
                return c;
            }        
            rs.close();
        }
        return null;
    }
    
    public Factor getFactor(String name) throws SQLException { 
        synchronized (loadFactorByName) {
            loadFactorByName.setString(1, name);
            ResultSet rs = loadFactorByName.executeQuery();
            
            if(rs.next()) { 
                Factor c = new Factor(rs);
                rs.close();
                
                if(!factorIDs.containsKey(c.getDBID())) { 
                    factorIDs.put(c.getDBID(), c);
                    factorNames.put(c.getName(), c);
                }
                
                return c;
            } else { 
                rs.close();
            }
        }
        int id = insertFactor(name);
        return loadFactor(id);
    }
    
    public Cells loadCells(int dbid) throws SQLException { 
        if(cellsIDs.containsKey(dbid)) { return cellsIDs.get(dbid); }

        Cells c = null;
        synchronized(loadCells) {
            loadCells.setInt(1, dbid);
            ResultSet rs = loadCells.executeQuery();
            if(rs.next()) { 
                c = new Cells(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Cells DBID: " + dbid);
            }
        }        
        cellsIDs.put(dbid, c);
        cellsNames.put(c.getName(), c);
        return c;
    }

    public Condition loadCondition(int dbid) throws SQLException { 
        if(condIDs.containsKey(dbid)) {  return condIDs.get(dbid); }
        
        Condition c = null;
        synchronized(loadCond) {
            loadCond.setInt(1, dbid);
            ResultSet rs = loadCond.executeQuery();
            
            if(rs.next()) { 
                c = new Condition(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Condition DBID: " + dbid);
            }
        }
        
        condIDs.put(dbid, c);
        condNames.put(c.getName(), c);
        return c;
    }

    public Factor loadFactor(int dbid) throws SQLException { 
        if(factorIDs.containsKey(dbid)) { return factorIDs.get(dbid); }
        
        
        Factor c = null;
        synchronized(loadFactor) {
            loadFactor.setInt(1, dbid);
            ResultSet rs = loadFactor.executeQuery();
            
            if(rs.next()) { 
                c = new Factor(rs);
                rs.close();
            } else {
                rs.close();
                throw new IllegalArgumentException("Unknown Condition DBID: " + dbid);
            }
        }
        
        factorIDs.put(dbid, c);
        factorNames.put(c.getName(), c);
        return c;
    }

    public Collection<Cells> loadAllCells(Collection<Integer> dbids) throws SQLException {

        LinkedList<Cells> values = new LinkedList<Cells>();
        for(int dbid : dbids) { values.addLast(loadCells(dbid)); }
        return values;
    }

    public Collection<Cells> loadAllCells() throws SQLException {
        
        HashSet<Cells> values = new HashSet<Cells>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id from cells");

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(loadCells(id));
        }

        rs.close();
        s.close();
        return values;
    }	
	
    public Collection<Condition> loadAllConditions(Collection<Integer> dbids) throws SQLException {

        LinkedList<Condition> values = new LinkedList<Condition>();
        for(int dbid : dbids) { values.addLast(loadCondition(dbid)); }
        return values;
    }

    public Collection<Condition> loadAllConditions() throws SQLException { 
        HashSet<Condition> values = new HashSet<Condition>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id from conditions");

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(loadCondition(id));
        }

        rs.close();
        s.close();
        return values;
    }	
	
    public Collection<Factor> loadAllFactors(Collection<Integer> dbids) throws SQLException {

        LinkedList<Factor> values = new LinkedList<Factor>();
        for(int dbid : dbids) { values.addLast(loadFactor(dbid)); }
        return values;
    }

    public Collection<Factor> loadAllFactors() throws SQLException { 
        HashSet<Factor> values = new HashSet<Factor>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id from factors");

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(loadFactor(id));
        }

        rs.close();
        s.close();
        return values;
    }	
	
    private int insertCells(String n) throws SQLException {
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select cells_id.nextval from dual");
        if(!rs.next()) { throw new IllegalArgumentException("cells_id sequence doesn't exist."); }
        int id = rs.getInt(1);
        rs.close();
        
        s.executeUpdate("insert into cells (id, name) values (" + id + ", '" + n + "')");
        
        s.close();
        return id;
    }

    private int insertCondition(String n) throws SQLException {
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select condition_id.nextval from dual");
        if(!rs.next()) { throw new IllegalArgumentException("condition_id sequence doesn't exist."); }
        int id = rs.getInt(1);
        rs.close();
        
        s.executeUpdate("insert into conditions (id, name) values (" + id + ", '" + n + "')");
        
        s.close();
        return id;
    }

    private int insertFactor(String n) throws SQLException {
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select factors_id.nextval from dual");
        if(!rs.next()) { throw new IllegalArgumentException("factors_id sequence doesn't exist."); }
        int id = rs.getInt(1);
        rs.close();
        
        s.executeUpdate("insert into factors (id, name) values (" + id + ", '" + n + "')");
        
        s.close();
        return id;
    }    
}
