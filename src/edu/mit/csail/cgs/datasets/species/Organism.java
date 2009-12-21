package edu.mit.csail.cgs.datasets.species;

import java.util.*;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import java.sql.*;

/**
 * <code>Organism</code> represents a species in our database
 */
public class Organism implements edu.mit.csail.cgs.utils.Closeable {

  private Map<String, Genome> genomes;

  private String species;

  private boolean isClosed;

  private int dbid;


  /**
   * 
   * @param species
   * @throws NotFoundException
   */
  public Organism(String species) throws NotFoundException {
    this.species = species;
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select id from species where name = '" + species + "'");
      
      if (rs.next()) {
          dbid = rs.getInt(1);
          rs.close();
          stmt.close();
      } else {
          DatabaseFactory.freeConnection(cxn);
          rs.close();
          stmt.close();
        throw new NotFoundException("Couldn't find " + species);
      } 

    }
    catch (SQLException ex) {
      ex.printStackTrace();
      throw new DatabaseException("Couldn't find " + species + ": " + ex.toString(), ex);
    }
    catch (UnknownRoleException ex) {
      ex.printStackTrace();
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
    if (genomes == null) {
      genomes = new HashMap<String, Genome>();
    }
    isClosed = false;
  }


  /**
   * 
   * @param speciesID
   * @throws NotFoundException
   */
  public Organism(int speciesID) throws NotFoundException {
    dbid = speciesID;
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select name from species where id=" + dbid);
      if (rs.next()) {
        species = rs.getString(1);
        rs.close();
        stmt.close();
      }
      else {
        rs.close();
        stmt.close();
        DatabaseFactory.freeConnection(cxn);
        throw new NotFoundException("Couldn't find " + dbid);
      }
    }
    catch (SQLException ex) {
      throw new DatabaseException("Couldn't find " + dbid + ": " + ex.toString(), ex);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
    if (genomes == null) {
      genomes = new HashMap<String, Genome>();
    }
    isClosed = false;
  }


  /* (non-Javadoc)
   * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
   */
  public boolean isClosed() {
    return isClosed;
  }


  /* (non-Javadoc)
   * @see edu.mit.csail.cgs.utils.Closeable#close()
   */
  public void close() {
    isClosed = true;
  }


  /**
   * @return
   */
  public String getName() {
    return species;
  }


  /**
   * @return
   */
  public int getDBID() {
    return dbid;
  }


  /**
   * @param version
   * @throws SQLException
   */
  public void insertGenome(String version) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement s = cxn.createStatement();
      String nextIdString = Sequence.getInsertSQL(cxn, "genome_id");
      String insertSQL = String.format("insert into genome(id, species, version) values (%s, %d, '%s')", nextIdString,
          dbid, version);
      s.executeUpdate(insertSQL);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    catch (SQLException se) {
      throw se;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }


  /**
   * Returns a Genome object for the specified build/version of this species.
   * @param version
   * @return
   * @throws NotFoundException
   */
  public Genome getGenome(String version) throws NotFoundException {
    String k = getName() + version;
    if (genomes.containsKey(k)) {
      return genomes.get(k);
    }
    else {
      Genome g = new Genome(dbid, version);
      genomes.put(k, g);
      return g;
    }
  }


  /**
   * Returns all of the versions/builds for this species.
   * @return
   */
  public Collection<String> getGenomeNames() {
    LinkedList<String> lst = new LinkedList<String>();
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select version from genome where species=" + dbid + " order by id");
      while (rs.next()) {
        lst.addLast(rs.getString(1));
      }
      rs.close();
      stmt.close();
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    catch (SQLException se) {
      se.printStackTrace(System.err);
    }
    DatabaseFactory.freeConnection(cxn);
    return lst;
  }


  /**
   * @return
   */
  public Collection<Integer> getGenomeIDs() {
    LinkedList<Integer> lst = new LinkedList<Integer>();
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select id from genome where species=" + dbid);
      while (rs.next()) {
        lst.addLast(rs.getInt(1));
      }
      rs.close();
      stmt.close();
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    catch (SQLException se) {
      se.printStackTrace(System.err);
    }
    DatabaseFactory.freeConnection(cxn);
    return lst;
  }

  static Map<String, Organism> organisms;

  static Map<Integer, Organism> organismsids;


  /**
   * @param genomeName
   * @return
   * @throws NotFoundException
   */
  public static Genome findGenome(String genomeName) throws NotFoundException {
    try {
      java.sql.Connection cxn = DatabaseFactory.getConnection("core");
      Statement s = cxn.createStatement();
      ResultSet rs = s.executeQuery("select species from genome where version='" + genomeName + "'");
      Genome g = null;

      if (rs.next()) {
        int orgID = rs.getInt(1);
        Organism org = getOrganism(orgID);
        try {
          g = org.getGenome(genomeName);
        }
        catch (NotFoundException nfe) {
          g = null;
        }
      }

      rs.close();
      s.close();
      DatabaseFactory.freeConnection(cxn);

      if (g == null) {
        throw new NotFoundException("Couldn't find genome: " + genomeName);
      }

      return g;

    }
    catch (SQLException se) {
      throw new DatabaseException("SQLException: " + se.getMessage(), se);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
  }


  /**
   * @param gid
   * @return
   * @throws NotFoundException
   */
  public static Genome findGenome(int gid) throws NotFoundException {
    try {
      java.sql.Connection cxn = DatabaseFactory.getConnection("core");
      Statement s = cxn.createStatement();
      ResultSet rs = s.executeQuery("select version, species from genome where id=" + gid);
      Genome g = null;

      if (rs.next()) {
        String genomeName = rs.getString(1);
        int orgID = rs.getInt(2);
        Organism org = getOrganism(orgID);
        try {
          g = org.getGenome(genomeName);
        }
        catch (NotFoundException nfe) {
          g = null;
        }
      }

      rs.close();
      s.close();
      DatabaseFactory.freeConnection(cxn);

      if (g == null) {
        throw new NotFoundException("Couldn't find genome: " + gid);
      }

      return g;

    }
    catch (SQLException se) {
      throw new DatabaseException("SQLException: " + se.getMessage(), se);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
  }


  /**
   * @param species
   * @return
   * @throws NotFoundException
   */
  public static Organism getOrganism(String species) throws NotFoundException {
    if (organisms == null) {
      organisms = new HashMap<String, Organism>();
      organismsids = new HashMap<Integer, Organism>();
    }
    if (organisms.containsKey(species)) {
      return organisms.get(species);
    }
    else {
      Organism o = new Organism(species);
      organisms.put(species, o);
      organismsids.put(o.getDBID(), o);
      return o;
    }
  }


  /**
   * @param species
   * @return
   * @throws NotFoundException
   */
  public static Organism getOrganism(int species) throws NotFoundException {
    if (organisms == null) {
      organisms = new HashMap<String, Organism>();
      organismsids = new HashMap<Integer, Organism>();
    }
    if (organismsids.containsKey(species)) {
      return organismsids.get(species);
    }
    else {
      Organism o = new Organism(species);
      organisms.put(o.getName(), o);
      organismsids.put(o.getDBID(), o);
      return o;
    }
  }


  /**
   * @param species
   * @throws SQLException
   */
  public static void insertOrganism(String species) throws SQLException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement s = cxn.createStatement();
      String nextIdString = Sequence.getInsertSQL(cxn, "species_id");
      String sql = String.format("insert into species values (%s, '%s')", nextIdString, species);
      System.out.println(sql);
      s.executeUpdate(sql);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    catch (SQLException se) {
      throw se;
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }


  /**
   * @return
   */
  public static Collection<String> getOrganismNames() {
    LinkedList<String> lst = new LinkedList<String>();
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select name from species");
      while (rs.next()) {
        lst.addLast(rs.getString(1));
      }
      rs.close();
      stmt.close();
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    catch (SQLException se) {
      se.printStackTrace(System.err);
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
    return lst;
  }


  /**
   * @param species
   * @return
   * @throws NotFoundException
   */
  public static int speciesID(String species) throws NotFoundException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select id from species where name ='" + species + "'");
      if (rs.next()) {
        return rs.getInt(1);
      }
      else {
        throw new NotFoundException("Can't find species " + species + " in database");
      }
    }
    catch (SQLException ex) {
      throw new DatabaseException("Can't get species id", ex);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }


  /**
   * @param species
   * @param version
   * @return
   * @throws NotFoundException
   */
  public static int genomeID(String species, String version) throws NotFoundException {
    return genomeID(speciesID(species), version);
  }


  /**
   * @param species
   * @param version
   * @return
   * @throws NotFoundException
   */
  public static int genomeID(int species, String version) throws NotFoundException {
    java.sql.Connection cxn = null;
    try {
      cxn = DatabaseFactory.getConnection("core");
      Statement stmt = cxn.createStatement();
      ResultSet rs = stmt.executeQuery("select id from genome where species = " + species + " and version = '"
          + version + "'");
      if (rs.next()) {
        return rs.getInt(1);
      }
      else {
        throw new NotFoundException("Can't find genome " + species + ", " + version);
      }
    }
    catch (SQLException ex) {
      throw new DatabaseException("Can't get genome id", ex);
    }
    catch (UnknownRoleException ex) {
      throw new DatabaseException("Couldn't connect with role core", ex);
    }
    finally {
      if (cxn != null) {
        DatabaseFactory.freeConnection(cxn);
      }
    }
  }
}
