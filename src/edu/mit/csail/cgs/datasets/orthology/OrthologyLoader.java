/*
 * Created on Oct 19, 2005
 */
package edu.mit.csail.cgs.datasets.orthology;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class OrthologyLoader implements edu.mit.csail.cgs.utils.Closeable {
    
    public static void main(String[] args) { 
        try {
            OrthologyLoader loader = new OrthologyLoader();
            for(OrthologyMapping m : loader.loadAllMappings()) { 
                System.out.println(m.getName() + "," + m.getVersion());
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public static String dbRole;
    
    static { 
        dbRole = "annotations";
    }
    
    private java.sql.Connection cxn;
    private PreparedStatement loadMappingID, loadMappingNV, loadMappingAll;
    private PreparedStatement loadPairID, loadPairFNG, loadPairSNG, loadPairMapping, loadPairMappingFNG, loadPairMappingSNG;
    
    private Map<Integer,OrthologyMapping> idMap;
    private Map<String,Map<String,OrthologyMapping>> nameVersionMap;

    public OrthologyLoader() throws SQLException {
        cxn = DatabaseFactory.getConnection(dbRole);
        loadMappingID = OrthologyMapping.prepareLoadByID(cxn);
        loadMappingNV = OrthologyMapping.prepareLoadByNameVersion(cxn);
        loadPairID = OrthologyPair.prepareLoadByID(cxn);
        loadPairFNG = OrthologyPair.prepareLoadByFirstNameGenome(cxn);
        loadPairSNG = OrthologyPair.prepareLoadBySecondNameGenome(cxn);
        loadPairMappingFNG = OrthologyPair.prepareLoadByFirstNameGenomeWithMapping(cxn);
        loadPairMappingSNG = OrthologyPair.prepareLoadBySecondNameGenomeWithMapping(cxn);
        loadPairMapping = OrthologyPair.prepareLoadByMapping(cxn);
        loadMappingAll = OrthologyMapping.prepareLoadAllStatement(cxn);

        idMap = new HashMap<Integer,OrthologyMapping>();
        nameVersionMap = new HashMap<String,Map<String,OrthologyMapping>>();
    }
    
    public Collection<OrthologyMapping> loadAllMappings() {
        try {
            LinkedList<OrthologyMapping> mappings = new LinkedList<OrthologyMapping>();
            synchronized (loadMappingAll) {
                ResultSet rs = loadMappingAll.executeQuery();
                
                while(rs.next()) { 
                    OrthologyMapping om = new OrthologyMapping(rs);
                    mappings.addLast(om);
                }
                rs.close();
            }
           
            return mappings;
        } catch (SQLException e) {
            throw new DatabaseException("Error: " + e.getMessage(), e);
        }
    }
    
    public OrthologyMapping loadMapping(int mappingID) throws NotFoundException { 
        if(idMap.containsKey(mappingID)) { return idMap.get(mappingID); }
        try { 
            boolean found = false;
            OrthologyMapping mapping = null;
            synchronized (loadMappingID) {
                loadMappingID.setInt(1, mappingID);
                ResultSet rs = loadMappingID.executeQuery();
                if(found = rs.next()) { 
                    mapping = new OrthologyMapping(rs);
                }
                rs.close();
            }
            if(!found) { throw new NotFoundException("Unknown Mapping ID: " + mappingID); }
            idMap.put(mappingID, mapping);
            if(!nameVersionMap.containsKey(mapping.getName())) { 
                nameVersionMap.put(mapping.getName(), new HashMap<String,OrthologyMapping>()); 
            }
            nameVersionMap.get(mapping.getName()).put(mapping.getVersion(), mapping);
            return mapping;
        } catch(SQLException se) { 
            throw new DatabaseException("SQLException occurred: " + se.getMessage(), se);
        }
    }
    
    public OrthologyMapping loadMapping(String name, String version) throws NotFoundException { 
        if(nameVersionMap.containsKey(name)) { 
            if(nameVersionMap.get(name).containsKey(version)) { 
                return nameVersionMap.get(name).get(version);
            }
        }
        
        try { 
            boolean found = false;
            OrthologyMapping mapping = null;
            synchronized (loadMappingNV) {
                loadMappingNV.setString(1, name);
                loadMappingNV.setString(2, version);
                ResultSet rs = loadMappingNV.executeQuery();
                if(found = rs.next()) { 
                    mapping = new OrthologyMapping(rs);
                }
                rs.close();
            }
            if(!found) { throw new NotFoundException("Unknown Mapping: " + name + "/" + version); }
            
            idMap.put(mapping.getDBID(), mapping);
            if(!nameVersionMap.containsKey(mapping.getName())) { 
                nameVersionMap.put(mapping.getName(), new HashMap<String,OrthologyMapping>()); 
            }
            nameVersionMap.get(mapping.getName()).put(mapping.getVersion(), mapping);
            return mapping;
        } catch(SQLException se) { 
            throw new DatabaseException("SQLException occurred: " + se.getMessage(), se);
        }
    }    
    
    public OrthologyPair getPair(int orthID) throws NotFoundException {
        try { 
            OrthologyPair p = null;
            synchronized (loadPairID) {
                loadPairID.setInt(1, orthID);
                ResultSet rs = loadPairID.executeQuery();
                if(rs.next()) { 
                    p = new OrthologyPair(rs, this);
                }
                rs.close();
            }
            if(p == null) { throw new NotFoundException("Unknown orthology pair: " + orthID); }
            return p;
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        }
    }
    
    public int countOrthologyPairs(OrthologyMapping mapping) throws NotFoundException { 
        try {
            Statement s = cxn.createStatement();
            int count = 0;
            ResultSet rs = s.executeQuery("select count(*) from orth_pair where mapping=" + mapping.getDBID());
            if(rs.next()) { 
                count = rs.getInt(1);
            }
            rs.close();
            s.close();
            return count;
        } catch (SQLException e) {
            throw new DatabaseException("Error: " + e.getMessage(), e);
        }
    }
    
    public TotalOrthologyMapping loadTotalOrthologyMapping(OrthologyMapping mapping) { 
        TotalOrthologyMapping tom = new TotalOrthologyMapping(mapping);

        try { 
            synchronized (loadPairMapping) {
                loadPairMapping.setInt(1, mapping.getDBID());
                ResultSet rs = loadPairMapping.executeQuery();
                while(rs.next()) { 
                    tom.addPair(new OrthologyPair(rs, this));
                }
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);            
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        
        return tom;
    }
    
    public Collection<OrthologyPair> getAllPairs(OrthologyMapping mapping) { 
        LinkedList<OrthologyPair> lst = new LinkedList<OrthologyPair>();
        try { 
            synchronized (loadPairMapping) {
                loadPairMapping.setInt(1, mapping.getDBID());
                ResultSet rs = loadPairMapping.executeQuery();
                while(rs.next()) { 
                    lst.addLast(new OrthologyPair(rs, this));
                }
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);            
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        return lst;
    }
    
    public Collection<OrthologyPair> getFirstNameGenomePairs(String name, Genome g) {
        LinkedList<OrthologyPair> lst = new LinkedList<OrthologyPair>();
        try {
            synchronized (loadPairFNG) {
                loadPairFNG.setString(1, name);
                loadPairFNG.setInt(2, g.getDBID());
                ResultSet rs = loadPairFNG.executeQuery();
                while(rs.next()) { 
                    lst.addLast(new OrthologyPair(rs, this));
                }
                
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        return lst;
    }

    public Collection<OrthologyPair> getSecondNameGenomePairs(String name, Genome g) {
        LinkedList<OrthologyPair> lst = new LinkedList<OrthologyPair>();
        try {
            synchronized (loadPairSNG) {
                loadPairSNG.setString(1, name);
                loadPairSNG.setInt(2, g.getDBID());
                ResultSet rs = loadPairSNG.executeQuery();
                while(rs.next()) { 
                    lst.addLast(new OrthologyPair(rs, this));
                }
                
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        return lst;
    }

    public Collection<OrthologyPair> getFirstNameGenomePairs(String name, Genome g, OrthologyMapping mapping) {
        LinkedList<OrthologyPair> lst = new LinkedList<OrthologyPair>();
        try {
            synchronized (loadPairMappingFNG) {
                loadPairMappingFNG.setString(1, name);
                loadPairMappingFNG.setInt(2, g.getDBID());
                loadPairMappingFNG.setInt(3, mapping.getDBID());
                ResultSet rs = loadPairMappingFNG.executeQuery();
                while(rs.next()) { 
                    lst.addLast(new OrthologyPair(rs, this));
                }
                
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        return lst;
    }

    public Collection<OrthologyPair> getSecondNameGenomePairs(String name, Genome g, OrthologyMapping mapping) {
        LinkedList<OrthologyPair> lst = new LinkedList<OrthologyPair>();
        try {
            synchronized (loadPairMappingSNG) {
                loadPairMappingSNG.setString(1, name);
                loadPairMappingSNG.setInt(2, g.getDBID());
                loadPairMappingSNG.setInt(3, mapping.getDBID());
                ResultSet rs = loadPairMappingSNG.executeQuery();
                while(rs.next()) { 
                    lst.addLast(new OrthologyPair(rs, this));
                }
                
                rs.close();
            }
        } catch(SQLException se) { 
            throw new DatabaseException("Error: " + se.getMessage(), se);
        } catch (NotFoundException e) {
            e.printStackTrace(System.err);
            throw new RuntimeException(e);
        }
        return lst;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        try {
            loadMappingNV.close(); loadMappingNV = null;
            loadMappingID.close(); loadMappingID = null;
            loadMappingAll.close(); loadMappingAll = null;

            loadPairID.close(); loadPairID = null;
            loadPairFNG.close(); loadPairFNG = null;
            loadPairSNG.close(); loadPairSNG = null;
            loadPairMappingFNG.close(); loadPairMappingFNG = null;
            loadPairMappingSNG.close(); loadPairMappingSNG = null;
            loadPairMapping.close(); loadPairMapping = null;

            DatabaseFactory.freeConnection(cxn);
        } catch(SQLException se) {
            throw new DatabaseException("Couldn't close connection: " + se.getMessage(), se);
        } finally { 
            cxn = null;
        }
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return cxn==null;
    }

}
