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
    PreparedStatement getSpeciesID, getGeneID, getGeneLinks, getGeneActions, getGeneName, getGenesLike, getGeneNameType;

    public StringDB() throws SQLException {
        cxn = DatabaseFactory.getConnection("string");
        getSpeciesID = cxn.prepareStatement("select id from species where compactname = ?");
        getGeneID = cxn.prepareStatement("select id from aliases where species = ? and alias = ?");
        getGenesLike = cxn.prepareStatement("select alias from aliases where species = ? and alias like '%' || ? || '%'");
        getGeneName = cxn.prepareStatement("select alias from aliases where id = ? and species = ?");
        getGeneLinks = cxn.prepareStatement("select geneB, score from proteinlinks where speciesA = ? and geneA = ? and speciesB = ?");
        getGeneActions = cxn.prepareStatement("select geneB, intmode, action, score from proteinactions where speciesA = ? and geneA = ? and speciesB = ?");
        getGeneNameType = cxn.prepareStatement("select alias from aliases where id = ? and species = ? and source like '%' || ? || '%'");
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
    public List<String> getGenesLike(int speciesID, String geneName) throws SQLException, NotFoundException {
        getGenesLike.setInt(1,speciesID);
        getGenesLike.setString(2,geneName);
        List<String> output = new ArrayList<String>();
        ResultSet rs = getGenesLike.executeQuery();
        while (rs.next()) {
            output.add(rs.getString(1));
        }
        rs.close();
        return output;
    }

    public String getGeneID(int speciesID, String geneName) throws SQLException, NotFoundException {
        getGeneID.setInt(1,speciesID);
        getGeneID.setString(2,geneName);
        ResultSet rs = getGeneID.executeQuery();
        if (rs.next()) {
            String result = rs.getString(1);
            rs.close();
            return result;
        } else {
            rs.close();
            throw new NotFoundException("No gene with name " + geneName + " and species " + speciesID);
        }
    }
    public List<String> getGeneName(int speciesID, String geneID) throws SQLException {
        getGeneName.setString(1,geneID);
        getGeneName.setInt(2,speciesID);
        ResultSet rs = getGeneName.executeQuery();
        List<String> output = new ArrayList<String>();
        while (rs.next()) {
            output.add(rs.getString(1));
        }
        rs.close();
        return output;
    }
    public List<String> getGeneName(int speciesID, String geneID, String type) throws SQLException {
        getGeneNameType.setString(1,geneID);
        getGeneNameType.setInt(2,speciesID);
        getGeneNameType.setString(3,type);
        ResultSet rs = getGeneNameType.executeQuery();
        List<String> output = new ArrayList<String>();
        while (rs.next()) {
            output.add(rs.getString(1));
        }
        rs.close();
        return output;        
    }
    public List<Link> getGeneLinks(int speciesID, String geneID) throws SQLException {
        ArrayList<Link> output = new ArrayList<Link>();
        getGeneLinks.setInt(1, speciesID);
        getGeneLinks.setString(2, geneID);
        getGeneLinks.setInt(3, speciesID);
        
        ResultSet rs = getGeneLinks.executeQuery();
        while (rs.next()) {
            output.add(new Link(speciesID, geneID, rs.getString(1), rs.getInt(2)));
        }
        rs.close();
        return output;
    }
    public List<Action> getGeneActions(int speciesID, String geneID) throws SQLException {
        ArrayList<Action> output = new ArrayList<Action>();
        getGeneActions.setInt(1, speciesID);
        getGeneActions.setString(2, geneID);
        getGeneActions.setInt(3, speciesID);
        ResultSet rs = getGeneActions.executeQuery();
        while (rs.next()) {
            output.add(new Action(speciesID, geneID, rs.getString(1), rs.getString(2), rs.getString(3), rs.getInt(4)));
        }
        rs.close();
        return output;
    }

    public void close() {
        if (cxn != null) {
            try {
                getSpeciesID.close();
                getGeneID.close();
                getGenesLike.close();
                getGeneName.close();
                getGeneNameType.close();
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