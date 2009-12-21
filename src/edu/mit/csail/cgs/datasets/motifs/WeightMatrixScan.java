package edu.mit.csail.cgs.datasets.motifs;
import java.util.*;
import java.sql.*;
import java.text.DecimalFormat;
import java.text.FieldPosition;

import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.NotFoundException;

public class WeightMatrixScan {
    public int scandbid;
    public boolean hasscandbid;
    public String scanname;
    public float cutoff;
    public WeightMatrix matrix;

    public static WeightMatrixScan getWeightMatrixScan(int dbid) throws NotFoundException {
        WeightMatrixScan scan = null;
        try {
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select weightmatrix,name,cutoff from weightmatrixscan where " +
                                                        " id = ?");
            ps.setInt(1,dbid);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                scan = new WeightMatrixScan();
                scan.scanname = rs.getString(2);
                scan.cutoff = rs.getFloat(3);
                scan.matrix = WeightMatrix.getWeightMatrix(rs.getInt(1));
                scan.scandbid = dbid;
                scan.hasscandbid = true;
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
            } else {
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
                throw new NotFoundException("Can't find WMSID " + dbid);
            }
        } catch (SQLException ex) {
            throw new NotFoundException("Can't find WMS " + dbid,ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Can't connect to annotations datasource",ex);
        }
        return scan;
    }

    public static Collection<WeightMatrixScan> getScansForMatrix (int matrixid) {
        java.sql.Connection cxn = null;
        ArrayList<WeightMatrixScan> results = new ArrayList<WeightMatrixScan>();
        try {
            cxn = DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select id from weightmatrixscan where weightmatrix = ?");
            ps.setInt(1,matrixid);
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                results.add(getWeightMatrixScan(rs.getInt(1)));
            }
            rs.close();
            ps.close();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }
        return results;
    }

    public static WeightMatrixScan getScanForMatrix(int matrixid, String scanname) {
        java.sql.Connection cxn = null;
        WeightMatrixScan results = null;
        try {
            cxn = DatabaseFactory.getConnection("annotations");
            PreparedStatement ps = cxn.prepareStatement("select id from weightmatrixscan where weightmatrix = ? and name = ?");
            ps.setInt(1,matrixid);
            ps.setString(2,scanname);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                results = getWeightMatrixScan(rs.getInt(1));
            }
            rs.close();
            ps.close();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }
        return results;
    }

    public static Collection<WeightMatrixScan> getScansForSpecies (int speciesid, String name, String version, String type, String scanname) {
        java.sql.Connection cxn = null;
        ArrayList<WeightMatrixScan> results = new ArrayList<WeightMatrixScan>();
        try {
            cxn = DatabaseFactory.getConnection("annotations");
            String sql = "select wms.id from weightmatrix wm, weightmatrixscan wms " +
                "where wm.id = wms.weightmatrix and wm.species = ?";
            if (name != null) { sql += " and wm.name = ? ";}
            if (version != null) {sql += " and wm.version = ? ";}
            if (scanname != null) {sql += " and wms.name = ? ";}
            if (type != null) {sql += " and wm.type = ? ";}
            sql += " order by wm.name, wm.version, wm.type";
            System.err.println("SQL " + sql);
            PreparedStatement ps = cxn.prepareStatement(sql);
            ps.setInt(1,speciesid);
            int count = 2;
            if (name != null) {ps.setString(count++,name);}
            if (version != null) { ps.setString(count++,version);}
            if (scanname != null) {ps.setString(count++,scanname);}
            if (type != null) {ps.setString(count++,type);}
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                results.add(getWeightMatrixScan(rs.getInt(1)));
            }
            rs.close();
            ps.close();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }
        return results;
    }

    public static Collection<WeightMatrixScan> getScansForGenome (int genomeid, String name, String version, String type, String scanname) {
        java.sql.Connection cxn = null;
        ArrayList<WeightMatrixScan> results = new ArrayList<WeightMatrixScan>();
        try {
            cxn = DatabaseFactory.getConnection("annotations");
            String sql = "select wms.id from weightmatrix wm, weightmatrixscan wms, wms_scanned_genomes wsg " +
                "where wm.id = wms.weightmatrix and wsg.genome = ? and wsg.scan = wms.id";
            if (name != null) { sql += " and wm.name = ? ";}
            if (version != null) {sql += " and wm.version = ? ";}
            if (scanname != null) {sql += " and wms.name = ? ";}
            if (type != null) {sql += " and wm.type = ? ";}
            sql += " order by wm.name, wm.version, wm.type";
            System.err.println("SQL " + sql);
            PreparedStatement ps = cxn.prepareStatement(sql);
            ps.setInt(1,genomeid);
            int count = 2;
            if (name != null) {ps.setString(count++,name);}
            if (version != null) { ps.setString(count++,version);}
            if (scanname != null) {ps.setString(count++,scanname);}
            if (type != null) {ps.setString(count++,type);}
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                results.add(getWeightMatrixScan(rs.getInt(1)));
            }
            rs.close();
            ps.close();
        } catch (NotFoundException ex) {
            ex.printStackTrace();
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }
        return results;
    }

    private static Collection<String> getField(int genomeid, String field) {
        java.sql.Connection cxn = null;
        ArrayList<String> results = new ArrayList<String>();
        try {
            cxn = DatabaseFactory.getConnection("annotations");
            String sql = "select unique(" + field + ") from weightmatrix wm, weightmatrixscan wms, wms_scanned_genomes wsg " +
                "where wm.id = wms.weightmatrix and wsg.genome = ? and wms.id = wsg.scan order by " + field;
            //            System.err.println("SQL IS " + sql);
            //            System.err.println("GENOME Is "+ genomeid);
            PreparedStatement ps = cxn.prepareStatement(sql);
            ps.setInt(1,genomeid);
            ResultSet rs = ps.executeQuery();
            while (rs.next()) {
                if (rs.getString(1) == null) {
                    //                    throw new NullPointerException("NULL getting " + field + " for " + speciesid);
                    System.err.println("  NOT ADDING " + rs.getString(1) + " for " + field + "," + genomeid);
                    continue;
                } 
                results.add(rs.getString(1));
            }
            rs.close();
            ps.close();
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
        } catch (SQLException ex) {
            ex.printStackTrace();
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }
        return results;
    }
    
    public static Collection<String> getNames(int genomeid) {return getField(genomeid,"wm.name");}
    public static Collection<String> getVersions(int genomeid) {return getField(genomeid,"wm.version");}
    public static Collection<String> getTypes(int genomeid) {return getField(genomeid,"wm.type");}
    public static Collection<String> getScanNames(int genomeid) {return getField(genomeid,"wms.name");}



    public String toString() {
        return "Scanned for (" + matrix.toString() + ") at " + cutoff;
    }
    public boolean equals(Object other) {
        if (other instanceof WeightMatrixScan) {
            WeightMatrixScan o = (WeightMatrixScan) other;
            return (o.matrix != null &&
                    matrix != null &&
                    o.matrix.equals(matrix) && 
                    scanname != null && 
                    o.scanname != null &&
                    scanname.equals(o.scanname) &&
                    cutoff == o.cutoff);
        } else {
            return false;
        }
    }
    public int hashCode() {
        return super.hashCode() + scanname.hashCode() + (int)(cutoff * 1900);
    }

}
