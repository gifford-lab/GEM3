package edu.mit.csail.cgs.datasets.motifs;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.NotFoundException;

public class WeightMatrixLoader implements edu.mit.csail.cgs.utils.Closeable {

    public WeightMatrixLoader() {}

    private Collection<String> queryWMTable(String field) {
        try {
            java.sql.Connection cxn = 
                DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select unique(" + field + ") from weightmatrix order by " + field);
            ResultSet rs = ps.executeQuery();
            ArrayList<String> out = new ArrayList<String>();
            while (rs.next()) {
                if (rs.getString(1) == null) {
                    System.err.println("Found a null string in " + field);
                } else {
                    out.add(rs.getString(1));
                }
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return out;
        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }
    }
    public Collection<String> getNames() {
        return queryWMTable("name");
    }
    public Collection<String> getVersions() {
        return queryWMTable("version");
    }

    public Collection<String> getTypes() {
        return queryWMTable("type");
    }

    public WeightMatrix query(int speciesid,
                              String name,
                              String version) {
        WeightMatrix wm = null;
        try {            
            java.sql.Connection cxn = 
                DatabaseFactory.getConnection("annotations");
            String query = "select wm.id from weightmatrix wm where wm.species = ? and " +
                    " wm.name = ? and wm.version = ?";
            PreparedStatement ps = cxn.prepareStatement(query);
            ps.setInt(1,speciesid);
            ps.setString(2,name);
            ps.setString(3,version);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                wm = WeightMatrix.getWeightMatrix(rs.getInt(1));
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);            
            if (wm != null) {
                return wm;
            } else {
                throw new NotFoundException ("Couldn't find " + name + "," + version + " in species " + speciesid);
            }

        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }
    
    }
    
    public Collection<WeightMatrix> loadMatrices(Organism species) throws SQLException { 
    	LinkedList<WeightMatrix> matrices = new LinkedList<WeightMatrix>();
    	
    	int speciesID = species.getDBID();
    	String wmquery = "select id, species, name, version, type from weightmatrix " +
    			"where species=?";
    	String colquery = "select position, letter, weight from weightmatrixcols where " +
    			"weightmatrix=?";
    	
        java.sql.Connection cxn = 
            DatabaseFactory.getConnection("annotations");
    	PreparedStatement wmStatement = cxn.prepareStatement(wmquery);
    	PreparedStatement colStatement = cxn.prepareStatement(colquery);
    	
    	wmStatement.setInt(1, speciesID);
    	ResultSet wmResults = wmStatement.executeQuery();
    	while(wmResults.next()) { 
    		int wmID = wmResults.getInt(1);
    		colStatement.setInt(1, wmID);
    		ResultSet colResults = colStatement.executeQuery();
    		
    		WeightMatrix matrix = new WeightMatrix(wmResults, colResults);
    		matrices.add(matrix);
    		
    		colResults.close();
    	}
    	wmResults.close();
    	
    	colStatement.close();
    	wmStatement.close();
    	DatabaseFactory.freeConnection(cxn);
    	
    	return matrices;
    }

    public Collection<WeightMatrix> query(String name,
                                          String version,
                                          String type) {
        if (name == null && version == null && type == null) {
            return WeightMatrix.getAllWeightMatrices();
        }
        try {
            java.sql.Connection cxn = 
                DatabaseFactory.getConnection("annotations");
            ArrayList<WeightMatrix> out = new ArrayList<WeightMatrix>();
            String query;
            if (name != null || version != null || type != null) {
                int index = 1;
                query = "select wm.id from weightmatrix wm where ";
                if (name != null) {
                    query += " wm.name = ? ";
                    index++;
                }
                if (version != null) {
                    if (index++ > 1) {
                        query += " and ";
                    }
                    query += " wm.version = ? ";
                }
                if (type != null) {
                    if (index++ > 1) {
                        query += " and ";
                    }
                    query += " wm.type = ? ";
                }
            } else {
                query = "select id from weightmatrix";
            }
            PreparedStatement ps = cxn.prepareStatement(query);
            if (name != null || version != null || type != null) {
                int index = 1;
                if (name != null) {
                    ps.setString(index++,name);
                }
                if (version != null) {
                    ps.setString(index++,version);
                }
                if (type != null) {
                    ps.setString(index++,type);
                }
            } 
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                WeightMatrix wm = WeightMatrix.getWeightMatrix(rs.getInt(1));
                if (wm != null) {
                    out.add(wm);
                }
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);

            return out;
        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (UnknownRoleException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        }

    }

    public void close() {}
    public boolean isClosed() {return true;}
}
