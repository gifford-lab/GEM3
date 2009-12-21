
package edu.mit.csail.cgs.datasets.chipchip;
import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public abstract class SQLGeneric implements GenericExperiment {
    private String chrom, tablename, name, version;
    /* some SQL* classes can either select all data in a window (chromosome, start, stop) or
       only the points for which the ratio exceeds some threshold.  lasttype
       indicates which type of query was performed last.   */
    private int windowstart, windowstop, genomeid, chromid, lasttype;
    private int[] positions;
    /* count is set by window() to be the number of items in the current window.
       exptid is the first experiment id that matches the supplied name, version, species.
       speciesid is the id of the species for the specified genomeid.
       sqlparambase is set by subclasses (in the constructor) to indiciate the index of the
         first parameter after the ones set in bindExptParams.   the index is 1-based*/
    protected int count, exptid, speciesid, sqlparambase;
    /* set by subclasses for use in window.  data retrieves the actual data that parseWindow expects.
       both should have parameters for the chromosome id, start, and stop after
       whatever parameters are set by bindExptParams. */
    protected String datasql;
    /* database connection.  protected so subclasses can use it */
    protected java.sql.Connection corecxn, chipcxn;
    private static Map<String,Integer> chromidmap;
    
    protected SQLGeneric(String tablename, String name, String version, int genomeid) throws NotFoundException {
        this.tablename = tablename;
        this.name = name;
        this.version = version;
        this.genomeid = genomeid;
        windowstart = -1;
        windowstop = -1;
        count = -1;
        lasttype = -1;
        try {
            corecxn = DatabaseFactory.getConnection("core");
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to database for role core",ex);
        }
        try {
            chipcxn = DatabaseFactory.getConnection("chipchip");
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to database for role chipchip",ex);
        }
        if (chromidmap == null) {
            chromidmap = new HashMap<String,Integer>();
        }
        /* default sqlparambase is 2 because there's normally one parameter set 
           by bindExptParams- the experiment id */
        sqlparambase = 2;
        try {
            /* get the species id */
            PreparedStatement stmt = corecxn.prepareStatement("select species from genome where id = ?");
            stmt.setInt(1,genomeid);
            ResultSet rs = stmt.executeQuery();
            if (rs.next()) {
                speciesid = rs.getInt(1);
            } else {
                throw new NotFoundException("Can't find genome with ID " + genomeid);
            }
            rs.close();
            stmt.close();
            /* get the experiment id */
            stmt = chipcxn.prepareStatement("select id from " + tablename + " where name = ? and species = ? and version = ?  and active = 1");
            stmt.setString(1,name);
            stmt.setInt(2,speciesid);
            stmt.setString(3,version);
            rs = stmt.executeQuery();
            if (rs.next()) {
                exptid = rs.getInt(1);
            } else {
				String msg = "Can't find exptid in " + tablename + " for (" + 
					name + "," + version + "," + speciesid + ")";
				throw new NotFoundException(msg);
                //throw new NotFoundException("Can'd find exptid for " + name + "," + version + " and " + speciesid + " in " + tablename);
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            System.err.println(tablename + "," + name + "," + version + "," + genomeid);
            throw new DatabaseException("Couldn't initialize SQLData " + ex.toString(),ex);
        }
    }
    /* turns chromosome names into ID for this species/genome.
       uses chromidmap to cache results */
    protected int getChromID() {return chromid;}
    protected int getChromID(String chrom) throws NotFoundException {
        String key = genomeid + "chr" + chrom;
        if (chromidmap.containsKey(key)) {
            return chromidmap.get(key);
        }
        try {
            PreparedStatement stmt = corecxn.prepareStatement("select id from chromosome where genome = ? and name = ?");
            stmt.setInt(1,genomeid);
            stmt.setString(2,chrom);
            ResultSet rs = stmt.executeQuery();
            if (rs.next()) {
                chromid = rs.getInt(1);
                rs.close();
                stmt.close();
                chromidmap.put(key,chromid);
                return chromid;
            } else {
                throw new NotFoundException("Can't find chromosome " + chrom + 
                                            "in genome " + genomeid);
            }
        } catch (SQLException ex) {
            throw new DatabaseException("Can't get chromid " + ex.toString(),ex);
        }        
    }
    protected void close() {
        DatabaseFactory.freeConnection(corecxn);
        DatabaseFactory.freeConnection(chipcxn);
    }
    public String getName() {return name;}
    public String getVersion() {return version;}
    protected int getGenomeID() {return genomeid;}
    protected void setPositions(int[] positions) {
        this.positions = positions;
    }
    protected void setWindowChrom(String c) {
        chrom = c;
    }
    protected void setWindowStart(int i) {
        windowstart = i;
    }
    protected void setWindowStop(int i) {
        windowstop = i;
    }
    public String getWindowChrom() {return chrom;}
    public int getWindowStart() {return windowstart;}
    public int getWindowStop() {return windowstop;}
    protected boolean sameWindow(String chrom, int start, int stop, int type) {
        if (chrom.equals(getWindowChrom()) &&
            start == getWindowStart() && 
            stop == getWindowStop() &&
            lasttype == type) {
            return true;
        }
        setWindowChrom(chrom);
        setWindowStart(start);
        setWindowStop(stop);
        lasttype = type;        
        return false;
    }
    /* sets the current data view to the described positions.

    This method uses datastmt and then calls parseWindow
    to interpret the results */
    public void window(String chrom, int start, int stop) throws NotFoundException {
        if (sameWindow(chrom,start,stop,1)) {return;}
        try {
            int chromid = getChromID(chrom);
            PreparedStatement datastmt = chipcxn.prepareStatement(datasql);
            datastmt.setFetchSize(1000);
            bindExptParams(datastmt);
            datastmt.setInt(sqlparambase + 0,chromid);
            datastmt.setInt(sqlparambase + 1,start);
            datastmt.setInt(sqlparambase + 2,stop);
            try {
                parseWindow(datastmt.executeQuery());
            } catch (ArrayIndexOutOfBoundsException ex) {
                System.err.println("Count is " + count);
                System.err.println("datasql is " + datasql);
                ex.printStackTrace();
                throw ex;
            }
            datastmt.close();
        } catch (SQLException ex) {
            System.err.println("Error in " + this + " trying to window");
            System.err.println("datastmt : " + datasql);
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    /* call PreparedStatment.set* to bind the parameters necessary to select
       the experiment.  In general, this sets the single parameter for the 
       experiment id.
    */
    protected void bindExptParams(PreparedStatement ps) throws SQLException {
        ps.setInt(1,exptid);
    }
    protected abstract void parseWindow(ResultSet rs) throws SQLException;
    
    public int getCount() {
        return count;
    }
    /* convert a base to the nearest index.  returns -1
       if the base isn't in this window */
    public int baseToIndex(int b) {
        if ((b > windowstop) || (b < windowstart)) {
            return -1;
        }
        try {
            int i = Arrays.binarySearch(positions,b);
            if (i < 0) {
                return -1*(i+1);
            } else {
                return i;
            }
        } catch (NullPointerException ex) {
            System.err.println("positions is null!!! " + positions);
            System.err.println("Name is " + getName() + " and version is " + getVersion());
            System.err.println("Window is " + getWindowChrom() +":" + 
                               getWindowStart() + "-" + getWindowStop());
            throw new RuntimeException("Null positions[]",ex);
        }
    }
    /* convert an index to a chromosomal position in this window*/
    public int getPos(int i) {
        return positions[i];
    }

    public String toString() {
        return getName()+"("+getVersion()+")["+getWindowChrom() +":" + 
            getWindowStart() + "-" + getWindowStop()+"]";
    }

    public double getMax(String chrom, int start, int stop) throws NotFoundException {
        window(chrom,start,stop);
        double max = (-1*Double.MAX_VALUE);
        for (int i = 0; i < getCount();i++) {
            for (int j = 0; j < getReplicates(i); j++) {
                if (getValue(i,j) > max) {max = getValue(i,j);}
            }
        }
        return max;
    }
    public double getMin(String chrom, int start, int stop) throws NotFoundException {
        window(chrom,start,stop);
        double min = Double.MAX_VALUE;
        for (int i = 0; i < getCount();i++) {
            for (int j = 0; j < getReplicates(i); j++) {
                if (getValue(i,j) < min) {min = getValue(i,j);}
            }
        }
        return min;
    }

}
