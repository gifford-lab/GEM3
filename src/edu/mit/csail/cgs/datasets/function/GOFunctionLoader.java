package edu.mit.csail.cgs.datasets.function;

/**
 * Loads functional annotations from a GO database.
 * Specify the database role (eg, go_200904) when
 * creating the loader.  Treats each species in the set of
 * annotations as a Version.
 *
 * FunctionVersion = species
 * Category = term
 * Assignment = association
 * Object = gene_product
 */

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.Closeable;

public class GOFunctionLoader implements FunctionLoader, Closeable {

    private java.sql.Connection cxn;
    private int is_a_dbid;

    public static String getDefaultDBName() {return "go_20100320";}

    public GOFunctionLoader(String goDBName) throws UnknownRoleException, SQLException {
        cxn = DatabaseFactory.getConnection(goDBName);
        Statement stmt = cxn.createStatement();
        ResultSet rs = stmt.executeQuery("select id from term where term_type = 'relationship' and name = 'is_a'");
        rs.next();
        is_a_dbid = rs.getInt(1);
        rs.close();
        stmt.close();
    }

	public Collection<FunctionVersion> getAllVersions() throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select id, concat(genus, ' ', species) from species");
        ResultSet rs = stmt.executeQuery();
        HashSet<FunctionVersion> output = new HashSet<FunctionVersion>();
        while (rs.next()) {
            output.add(new FunctionVersion(rs));
        }
        rs.close();
        stmt.close();
        return output;
    }

    public FunctionVersion getVersion(String versionName)  throws SQLException {
        /* assume the genus has no spaces in it.  This is 
           true as of the 200904 load
        */
        int splitAt = versionName.indexOf(" ");
        String genus = versionName.substring(0,splitAt);
        String species = versionName.substring(splitAt+1);
        PreparedStatement stmt = cxn.prepareStatement("select id, concat(genus, ' ', species) from species where genus = ? and species = ?");
        stmt.setString(1,genus);
        stmt.setString(2,species);
        ResultSet rs = stmt.executeQuery();
        rs.next();
        FunctionVersion output = new FunctionVersion(rs);
        rs.close();
        stmt.close();
        return output;
    }
    public FunctionVersion getVersion(int dbid) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select id, concat(genus, ' ', species) from species where id = ? ");
        stmt.setInt(1,dbid);
        ResultSet rs = stmt.executeQuery();
        rs.next();
        FunctionVersion output = new FunctionVersion(rs);
        rs.close();
        stmt.close();
        return output;
    }
    public Category getCategory(int dbid) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select id, name from term where id = ? ");
        stmt.setInt(1,dbid);
        ResultSet rs = stmt.executeQuery();
        rs.next();
        Category output = new Category(rs,this);
        rs.close();
        stmt.close();
        return output;
    }

    public Category getCategory(FunctionVersion fv, String name) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select term.id, term.name from term where term.name = ?");
        stmt.setString(1,name);
        ResultSet rs = stmt.executeQuery();
        rs.next();
        Category output = new Category(fv, name, name, rs.getInt(1));
        rs.close();
        stmt.close();
        return output;
    }

	public Collection<Category> getCategories(FunctionVersion fv) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select term.id, term.name from term where term.id in (select term_id from association where gene_product_id in (select id from gene_product where species_id = ?))");
        stmt.setInt(1,fv.getID());
        HashSet<Category> output = new HashSet<Category>();
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            output.add(new Category(fv,
                                    rs.getString(2), 
                                    rs.getString(2),
                                    rs.getInt(1)));
        }
        rs.close();
        stmt.close();
        return output;
    }

    public Collection<Assignment> getAssignments(FunctionVersion version) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select gene_product.symbol, association.term_id from association, gene_product where gene_product.id = association.gene_product_id and gene_product.species_id = ?");
        stmt.setInt(1,version.getID());
        HashSet<Assignment> output = new HashSet<Assignment>();
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            output.add(new Assignment(rs,this));
        }
        rs.close();
        stmt.close();
        return output;
    }

    public Collection<Assignment> getAssignments(Category c) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select gene_product.symbol, association.term_id from association, gene_product where gene_product.id = association.gene_product_id and association.term_id = ?");
        stmt.setInt(1,c.getID());
        HashSet<Assignment> output = new HashSet<Assignment>();
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            output.add(new Assignment(rs,this));
        }
        rs.close();
        stmt.close();
        return output;        
    }

	public Collection<Assignment> getAssignments(String obj, FunctionVersion fv) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select gene_product.symbol, association.term_id from association, gene_product where gene_product.id = association.gene_product_id and gene_product.species_id = ? and gene_product.symbol = ?");
        stmt.setInt(1,fv.getID());
        stmt.setString(2,obj);
        HashSet<Assignment> output = new HashSet<Assignment>();
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            output.add(new Assignment(rs,this));
        }
        rs.close();
        stmt.close();
        return output;
    }

	public Collection<Assignment> getAssignments(Category c, FunctionVersion fv) throws SQLException {
        PreparedStatement stmt = cxn.prepareStatement("select gene_product.symbol from association, gene_product where gene_product.id = association.gene_product_id and gene_product.species_id = ? and association.term_id = ?");
        stmt.setInt(1,fv.getID());
        stmt.setInt(2,c.getID());
        HashSet<Assignment> output = new HashSet<Assignment>();
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            output.add(new Assignment(rs.getString(1),c));
        }
        rs.close();
        stmt.close();
        return output;
    }
    
    public Collection<Assignment> getAllAssignments(Category c) throws SQLException {
        return getAssignments(c);
    }
    
    public Collection<Category> getChildCategories(Category c) throws SQLException {
        int isa = -1, partof = -1;
        PreparedStatement stmt = cxn.prepareStatement("select id,name from term where is_relation = 1 and (name = 'is_a' or name = 'part_of')");
        ResultSet rs = stmt.executeQuery();
        while (rs.next()) {
            if (rs.getString(2).equals("is_a")) {
                isa = rs.getInt(1);
            } else if (rs.getString(2).equals("part_of")) {
                partof = rs.getInt(1);
            }
        }

        stmt = cxn.prepareStatement("select term2term.term2_id, term.name from term, term2term where term2term.term2_id = term.id and term2term.term1_id = ? and (term2term.relationship_type_id = " +
                                    isa + " or relationship_type_id = " + partof + ")");
        stmt.setInt(1,c.getID());
        HashSet<Category> output = new HashSet<Category>();
        rs = stmt.executeQuery();
        while (rs.next()) {
            /* pass null as the FunctionVersion because GO terms aren't necessarily associated
               with a single species or such.  */               
            output.add(new Category(null, rs.getString(2),
                                    rs.getString(2),
                                    rs.getInt(1)));
                                    
        }
        rs.close();
        stmt.close();
        return output;       
    }

    

    public void close() {
        if (cxn != null) {
            DatabaseFactory.freeConnection(cxn);
            cxn = null;
        }
    }
    public boolean isClosed() {return cxn == null;}
}