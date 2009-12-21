/*
 * Created on Oct 20, 2005
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
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author tdanford
 */
public class OrthologyParser implements edu.mit.csail.cgs.utils.Closeable {
    
    private java.sql.Connection cxn;
    private PreparedStatement mappingInsert, pairInsert, mappingDelete, mappingPairDelete, pairDelete;
    private Map<String,Genome> genomeCache;

    public OrthologyParser(String user, String pword) throws SQLException, UnknownRoleException {
        cxn = DatabaseFactory.getConnection(OrthologyLoader.dbRole, user, pword);
        mappingInsert = OrthologyMapping.prepareInsertStatement(cxn);
        pairInsert = OrthologyPair.prepareInsertStatement(cxn);
        mappingDelete = OrthologyMapping.prepareDeleteStatement(cxn);
        mappingPairDelete = OrthologyMapping.preparePairDeleteStatement(cxn);
        pairDelete = OrthologyPair.prepareDeleteStatement(cxn);
        
        genomeCache = new HashMap<String,Genome>();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        try {
            mappingInsert.close();
            pairInsert.close();
            mappingDelete.close();
            mappingPairDelete.close();
            pairDelete.close();
            cxn.close();
            genomeCache.clear();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        cxn = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return cxn == null;
    }
    
    public TotalOrthologyMapping parseTotalOrthologyMapping(File f) throws IOException { 
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        String name = br.readLine();
        String version = br.readLine();
        if(name == null || version == null) { 
            br.close();
            throw new IOException("Name, Version weren't present at start of file.");
        }
        OrthologyMapping m = new OrthologyMapping(name, version);
        TotalOrthologyMapping tom = new TotalOrthologyMapping(m);
        int count = 0;
        
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0) { 
                String[] array = line.split("[\\s]+");
                String n1 = array[0], n2 = array[2];
                String g1 = array[1], g2 = array[3];
                try { 
                    if(!genomeCache.containsKey(g1)) { 
                        Genome gc = Organism.findGenome(g1);
                        genomeCache.put(g1, gc);
                    }
                    if(!genomeCache.containsKey(g2)) { 
                        Genome gc = Organism.findGenome(g2);
                        genomeCache.put(g2, gc);
                    }
                } catch(NotFoundException nfe) { 
                    br.close();
                    throw new IOException("Unknown Genome(s): " + g1 + "/" + g2);
                }
                
                Genome genome1 = genomeCache.get(g1);
                Genome genome2 = genomeCache.get(g2);
                OrthologyPair op = new OrthologyPair(m, n1, genome1, n2, genome2);
                tom.addPair(op);
                count += 1;
            }
        }
        
        br.close();
        System.out.println("# Orthology Pairs: " + count);
        return tom;
    }
    
    public void deleteOrthologyMapping(OrthologyMapping mapping) throws SQLException {
        if(mapping.getDBID() == -1) { 
            throw new IllegalArgumentException("Can't delete non-DB OrthologyMapping from DB"); 
        }
        try {
            cxn.setAutoCommit(false);
            mappingPairDelete.setInt(1, mapping.getDBID());
            mappingPairDelete.executeUpdate();
            mapping.deleteFromDB(mappingDelete);
            cxn.commit();
        } catch(SQLException se) { 
            throw se;
        } finally { 
            cxn.setAutoCommit(true);
        }
    }
    
    public void insertTotalOrthologyMapping(TotalOrthologyMapping tom) throws SQLException { 
        cxn.setAutoCommit(false);
        int nextMappingID = findMaxMappingID() + 1;
        int nextPairID = findMaxPairID() + 1;
        tom.insertIntoDB(mappingInsert, pairInsert, nextMappingID, nextPairID);
        cxn.commit();
        cxn.setAutoCommit(true);
    }

    public int findMaxMappingID() throws SQLException { 
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select max(id) from orth_mapping");
        int maxID = -1;
        if(rs.next()) { 
            maxID = rs.getInt(1);
        }
        rs.close();
        s.close();
        return maxID;
    }

    public int findMaxPairID() throws SQLException { 
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select max(id) from orth_pair");
        int maxID = -1;
        if(rs.next()) { 
            maxID = rs.getInt(1);
        }
        rs.close();
        s.close();
        return maxID;
    }
}
