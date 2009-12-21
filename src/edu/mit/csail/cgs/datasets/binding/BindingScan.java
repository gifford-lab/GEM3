package edu.mit.csail.cgs.datasets.binding;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.chipchip.NameVersion;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/**
 * @author Timothy Danford
 *
 * The objects of the BindingScan class represent rows from the "bindingscan" table.
 * 
 * The one wrinkle is: this isn't *exactly* correct.  BindingScan objects actually also carry
 * a Genome object with them; that is, they are rows from the bindingscan table associated with
 * a particular genome in the DB.  
 * 
 * It is *not* required that the Genome they are associated with also show up in the 
 * bindingscanToGenome table, although it usually will.  All that is required is that the 
 * associated experiment/bayes/MSP be mapped to the given genome.
 */

public class BindingScan {
	
    public static final int AGILENT_TYPE = 0;
    public static final int BAYES_TYPE = 1;
    public static final int MSP_TYPE = 2;

    private int dbid;
    private Genome genome;
    private String type, version;

    /**
     * There are two constructors for BindingScan.  This constructor takes most of the fields
     * as arguments, explicitly, and creates an out-of-database BindingScan object.
     * This object should then be inserted into the database via the BindingScanLoader.insertScan()
     * method, at which point it will be assigned a valid DBID.
     * 
     * @param g
     * @param t
     * @param v
     * @param et
     * @param exptArray
     * @param prms
     * @param rs
     */
    public BindingScan(Genome g, String t, String v) { 
        dbid = -1;
        genome = g;
        type = t;
        version = v;        
    }
    
    /**
     * This constructor, which should be only package-accessible, is called by
     * BindingScanLoader as it loads BindingScan objects from the database.
     * 
     * @param g
     * @param rs
     * @param paramStatement
     * @param regionStatement
     * @throws SQLException
     */
    BindingScan(Genome g, ResultSet rs) throws SQLException {
		
        genome = g;
		
        dbid = rs.getInt(1);
        type = rs.getString(2);
        version = rs.getString(3);        
    }

    /*
     * These are the basic accessors for the BindingScan object.
     */
	
    public int getDBID() { return dbid; }
    public String getType() { return type; }
    public String getVersion() { return version; }
    public Genome getGenome() { return genome; }
    
    public String toString() { return version + " (" + type + ")"; }
	
    public int hashCode() { 
        int code = 17;
        code += genome.hashCode(); code *= 37;
        code += version.hashCode(); code *= 37;
        code += type.hashCode(); code *= 37;
        return code;
    }
	
    public boolean equals(Object o) { 
        if(!(o instanceof BindingScan)) { return false; }
        BindingScan bs = (BindingScan)o;
        if(!genome.equals(bs.genome)) { return false; }
        if(!version.equals(bs.version)) { return false; }
        if(!type.equals(bs.type)) { return false; }
        return true;
    }
	
    /**
     * insertIntoDB is called by the BindingScanLoader object, when it wants to 
     * put this BindingScan into the database.  Users who want to insert a BindingScan
     * object should call the BindingScanLoader.insertScan() method.
     * 
     * @param c
     * @param genomeMap
     * @param ps
     * @param paramps
     * @param regionps
     * @throws SQLException
     */
    void insertIntoDB(java.sql.Connection c,
                      BindingScanGenomeMap genomeMap,
                      PreparedStatement ps)
        throws SQLException {
        
        if(dbid != -1) { throw new IllegalArgumentException(); }
        
        ps.setString(1, type);
        ps.setString(2, version);
		
        ps.executeUpdate();
        
        Statement s = c.createStatement();
        ResultSet rs = s.executeQuery("select bindingscan_id.currval from dual");
        if(rs.next()) { dbid = rs.getInt(1); }
        rs.close();
        s.close();
        
        genomeMap.insertNewMapping(c, dbid, genome.getDBID());
    }
	
    /*
     * The following methods prepare particular Statement objects for rapid insertion 
     * and deletion of BindingScan objects and their associated data into the database.
     * They are primarily used by the BindingScanLoader class. 
     */
	
    public static PreparedStatement prepareInsertStatement(java.sql.Connection c) throws SQLException { 
        String str = "insert into bindingscan (id, type, version) values " +
            "(bindingscan_id.nextval, ?, ?)";
        return c.prepareStatement(str);
    }
	
    public static PreparedStatement prepareInsertParamStatement(java.sql.Connection c) 
        throws SQLException { 
        String str = "insert into bindingscanparam (scan, key, value) values (?, ?, ?)";
        return c.prepareStatement(str);
    }
    
    public static PreparedStatement prepareInsertRegionStatement(java.sql.Connection c) throws SQLException { 
        String str = "insert into bindingscanregion (scan, chromosome, startpos, stoppos) values (?, ?, ?, ?)";
        return c.prepareStatement(str);
    }
    
    public static PreparedStatement prepareInsertExptStatement(java.sql.Connection c) throws SQLException { 
        String str = "insert into bindingscanToExpt (scan, scanexpt, scantype) values (?, ?, ?)";
        return c.prepareStatement(str);
    }
    
    public static PreparedStatement prepareLoadByID(java.sql.Connection c) throws SQLException {
        String query = "select id, type, version from " +
            "bindingscan where id=?";
        return c.prepareStatement(query);
    }

    public static PreparedStatement prepareLoadByRegion(java.sql.Connection c) throws SQLException {
        String str = "select s.id, s.type, s.version from bindingscan s where " + 
            " s.id in (select r.scan from bindingscanregion r where r. chromosome=? and " +
            "((r.startpos <= ? and r.stoppos >= ?) or (r.startpos >= ? and r.startpos <= ?)))";
        return c.prepareStatement(str);
    }
    
    public static PreparedStatement prepareLoadByVersionType(java.sql.Connection c) throws SQLException { 
        String query = "select id, type, version from " +
            "bindingscan where version=? and type=?";
        return c.prepareStatement(query);
    }
    
    public static PreparedStatement prepareLoadByLikeVersionType(java.sql.Connection c) throws SQLException { 
        String query = "select id, type, version from " +
            "bindingscan where version like ? and type=?";
        return c.prepareStatement(query);
    }
    
    /*
     * Finally, we get a series of helper methods, used here and elsewhere.  Included 
     * is a method (getLocatorType) which will take an ExptLocator, and return the array 
     * of int's that correspond to that Locator (useful, if you're initializing a BindingScan
     * object for the first time).
     */
	
    public static int getLocatorType(ExptLocator loc) { 
        if(loc instanceof ChipChipLocator) { return AGILENT_TYPE; }
        if(loc instanceof BayesLocator) { return BAYES_TYPE; }
        if(loc instanceof MSPLocator) { return MSP_TYPE; }
        return -1;
    }
	
    public static int[] getLocatorIDs(ExptLocator loc, Connection c) throws SQLException {
        int[] array = null;
        if(loc instanceof ChipChipLocator) { array = getChipChipIDs((ChipChipLocator)loc, c); }
        if(loc instanceof BayesLocator) { array = getBayesIDs((BayesLocator)loc, c); }
        if(loc instanceof MSPLocator) { array = getMSPIDs((MSPLocator)loc, c); }		
        return array;
    }
	
    public static int[] getChipChipIDs(ChipChipLocator loc, Connection c) throws SQLException { 
        Statement s = c.createStatement();
        String name = loc.getNameVersion().name;
        String version = loc.getNameVersion().version;
        ResultSet rs = s.executeQuery("select id from experiment where name='" + name + "' and " +
                                      "version = '" + version + "'");
        Vector<Integer> values = new Vector<Integer>();
        while(rs.next()) { 
            values.add(rs.getInt(1));
        }
        s.close();
        int[] array = new int[values.size()];
        for(int i = 0; i < values.size(); i++) { 
            array[i] = values.get(i);
        }
        Arrays.sort(array);
        return array;
    }

    public static int[] getBayesIDs(BayesLocator loc, Connection c) throws SQLException { 
        Statement s = c.createStatement();
        String name = loc.getNameVersion().name;
        String version = loc.getNameVersion().version;
        ResultSet rs = s.executeQuery("select id from bayesanalysis where name='" + name + "' and " +
                                      "version = '" + version + "'");
        Vector<Integer> values = new Vector<Integer>();
        while(rs.next()) { 
            values.add(rs.getInt(1));
        }
        s.close();
        int[] array = new int[values.size()];
        for(int i = 0; i < values.size(); i++) { 
            array[i] = values.get(i);
        }
        Arrays.sort(array);
        return array;
    }

    public static int[] getMSPIDs(MSPLocator loc, Connection c) throws SQLException { 
        Statement s = c.createStatement();
        String name = loc.getNameVersion().name;
        String version = loc.getNameVersion().version;
        ResultSet rs = s.executeQuery("select id from rosettaanalysis where name='" + name + "' and " +
                                      "version = '" + version + "'");
        Vector<Integer> values = new Vector<Integer>();
        while(rs.next()) { 
            values.add(rs.getInt(1));
        }
        s.close();
        int[] array = new int[values.size()];
        for(int i = 0; i < values.size(); i++) { 
            array[i] = values.get(i);
        }
        Arrays.sort(array);
        return array;
    }
}

/*

create sequence bindingscan_id;
create table bindingscan (
id number(10) unique not null,
version varchar2(200) not null,
type varchar2(100) not null
);

create table bindingscanToExpt (
scan constraint fk_bindingscantoexpt_scan references bindingscan(id) not null,
scanexpt varchar2(200) not null,
scantype number(10) not null
);

create table bindingscanToGenome (
scan constraint fk_bindingscantogenome_scan references bindingscan(id) not null,
genome constraint fk_bindingscantogenome_genome references genome(id) not null
);

create table bindingscanregion (
scan constraint fk_bindingscanregion_scan references bindingscan(id) not null,
chromosome constraint fk_bindingscan_chromosome references chromosome(id) not null,
startpos number(10) not null,
stoppos number(10) not null,
constraint bindingscanregion_pk primary key(scan,chromosome,startpos,stoppos)
) organization index compress 2;

create table bindingscanparam (
scan constraint fk_bsp_bindingscan references bindingscan(id) not null,
key varchar2(50),
value varchar2(50)
);

create table bindingevent (
scan constraint fk_bindingevent_scan references bindingscan(id) not null,
chromosome constraint fk_bindingevent_chromosome references chromosome(id) not null,
startpos number(10) not null,
stoppos number(10) not null,
eventsize binary_float not null,
eventconf binary_float not null,
constraint bindingevent_pk primary key (scan, chromosome, startpos, stoppos)
) organization index compress 2;

*/

