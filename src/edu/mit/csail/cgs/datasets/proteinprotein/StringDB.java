package edu.mit.csail.cgs.datasets.proteinprotein;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.DatabaseException;


/**
 * An interface to our copy of the the String protein-protein interaction database. 
 * http://string-db.com/
 *
 */

public class StringDB implements edu.mit.csail.cgs.utils.Closeable {

    private java.sql.Connection cxn;
    PreparedStatement getSpeciesID, getGeneID, getGeneLinks, getGeneActions, getGeneName;

    public StringDB() throws SQLException {
        cxn = DatabaseFactory.getConnection("string");
        getSpeciesID = cxn.prepareStatement("select id from species where compactname = ?");
        getGeneID = cxn.prepareStatement("select id from aliases where species = ? and alias = ?");
        getGeneName = cxn.prepareStatement("select alias from aliases where id = ?");
        getGeneLinks = cxn.prepareStatement("select b, score from proteinlinks where a = ?");
        getGeneActions = cxn.prepareStatement("select b, intmode, action, score from proteinactions where a = ?");       
    }
    public Integer getSpeciesID(Genome g) throws SQLException, NotFoundException {
        return getSpeciesID(g.getSpecies());
    }
    public Integer getSpeciesID(Organism o) throws SQLException, NotFoundException {
        return getSpeciesID(o.getName());
    }
    /** Returns the ID for the species */
    public int getSpeciesID(String species) throws SQLException, NotFoundException {
        getSpeciesID.setString(1,species);
        ResultSet rs = getSpeciesID.executeQuery();
        if (rs.next()) {
            int result = rs.getInt(1);
            rs.close();
            return result;
        } else {
            rs.close();
            throw new NotFoundException("No species in stringdb for " + species);
        }
    }

    public int getGeneID(String geneName, int speciesID) throws SQLException, NotFoundException {
        getGeneID.setInt(1,speciesID);
        getGeneID.setString(2,geneName);
        ResultSet rs = getGeneID.executeQuery();
        if (rs.next()) {
            int result = rs.getInt(1);
            rs.close();
            return result;
        } else {
            rs.close();
            throw new NotFoundException("No gene with name " + geneName + " and species " + speciesID);
        }
    }
    public String getGeneName(int geneID) throws SQLException, NotFoundException {
        getGeneName.setInt(1,geneID);
        ResultSet rs = getGeneName.executeQuery();
        if (rs.next()) {
            String name = rs.getString(1);
            rs.close();
            return name;
        } else {
            rs.close();
            throw new NotFoundException("No such gene with id " + geneID);
        }
    }
    public List<Link> getGeneLinks(int geneID) throws SQLException {
        ArrayList<Link> output = new ArrayList<Link>();
        getGeneLinks.setInt(1,geneID);
        ResultSet rs = getGeneLinks.executeQuery();
        while (rs.next()) {
            output.add(new Link(geneID, rs.getInt(1), rs.getInt(2)));
        }
        rs.close();
        return output;
    }
    public List<Action> getGeneActions(int geneID) throws SQLException {
        ArrayList<Action> output = new ArrayList<Action>();
        getGeneActions.setInt(1,geneID);
        ResultSet rs = getGeneActions.executeQuery();
        while (rs.next()) {
            output.add(new Action(geneID, rs.getInt(1), rs.getString(2), rs.getString(3), rs.getInt(4)));
        }
        rs.close();
        return output;
    }

    public void close() {
        if (cxn != null) {
            try {
                getSpeciesID.close();
                getGeneID.close();
                getGeneName.close();
                getGeneActions.close();
                getGeneLinks.close();
                cxn.close();
            } catch (SQLException e) {
                throw new DatabaseException(e.toString(), e);
            }
            cxn = null;
        }

    }
    public boolean isClosed() {return cxn == null;}
}