package edu.mit.csail.cgs.datasets.binding;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

public class BindingScanLoader implements edu.mit.csail.cgs.utils.Closeable {
    
    public static void main(String[] args) { 
        try {
            BindingScanLoader loader = new BindingScanLoader();
            
            if(args.length > 0 && args[0].equals("delete")) { 
            	if(args[1].equals("all")) {
            		LinkedList<Integer> dbids = new LinkedList<Integer>();
            		java.sql.Connection cxn = loader.getConnection();
            		Statement s= cxn.createStatement();
            		ResultSet rs = s.executeQuery("select id from bindingscan");
            		while(rs.next()) { dbids.addLast(rs.getInt(1)); }
            		rs.close();
            		s.close();
            		
            		for(int dbid : dbids) { 
            			loader.deleteScan(dbid);
            		}
            	} else { 
            		int dbid = Integer.parseInt(args[1]);
            		loader.deleteScan(dbid);
            	}
            }
            
            if(args.length == 0 || args[0].equals("list")) {
                String gname = args.length > 1 ? args[1] : "sacCer1";
                Genome g = Organism.findGenome(gname);
                Collection<BindingScan> scans = loader.loadScans(g);
                for(BindingScan scan : scans) {
                    System.out.println(scan.getVersion() + " // " + scan.getType());
                }
            }
            
            loader.close();

        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public static final String role = "chipchip";

	private java.sql.Connection c;
    private BindingScanGenomeMap genomeMap;
	private PreparedStatement loadScanByRegion, loadScanByID, loadScansByVersionType, loadParams, loadRegions, loadExpts;
	private PreparedStatement loadEventByScan, loadEventByScanRegion;
	private PreparedStatement insertScan, insertParam, insertRegion, insertEvent, insertExpt;
	
	public BindingScanLoader() throws SQLException, UnknownRoleException { 
		c = DatabaseFactory.getConnection(role);
        init();
	}
    
    private void init() throws SQLException, UnknownRoleException { 
        genomeMap = new BindingScanGenomeMap(c);
        
        loadScanByRegion = BindingScan.prepareLoadByRegion(c);
        loadScansByVersionType = BindingScan.prepareLoadByVersionType(c);
        loadScanByID = BindingScan.prepareLoadByID(c);
        loadParams = loadParamsByScan(c);
        loadRegions = loadRegionsByScan(c);
        loadExpts = loadExptsByScan(c);
        loadEventByScan = prepareLoadEventByScan(c);
        loadEventByScanRegion = prepareLoadEventByScanRegion(c);
        
        insertScan = BindingScan.prepareInsertStatement(c);
        insertParam = BindingScan.prepareInsertParamStatement(c);
        insertRegion = BindingScan.prepareInsertRegionStatement(c);
        insertExpt = BindingScan.prepareInsertExptStatement(c);
        insertEvent = prepareInsertEvent(c);
    }
    
    public java.sql.Connection getConnection() { return c; }

	public void close() { 
		try {
			loadScanByRegion.close(); loadScanByRegion = null;
            loadScansByVersionType.close(); loadScansByVersionType = null;
			loadScanByID.close(); loadScanByID = null;
			loadParams.close(); loadParams = null;
            loadRegions.close(); loadRegions = null;
            loadExpts.close(); loadExpts = null;
			loadEventByScan.close(); loadEventByScan = null;
            loadEventByScanRegion.close(); loadEventByScanRegion = null;
			
			insertScan.close(); insertScan = null;
			insertParam.close(); insertParam = null; 
            insertRegion.close(); insertRegion = null;
			insertEvent.close(); insertEvent = null;
            insertExpt.close(); insertExpt = null;
			
			DatabaseFactory.freeConnection(c);
		} catch (SQLException e) {
			e.printStackTrace();
		}

        c = null;
	}
	
	public boolean isClosed() { 
		return c == null;
	}
    
    public Collection<BindingScan> loadScans(Genome g) throws SQLException { 
        if(!genomeMap.containsGenome(g.getDBID())) { return new LinkedList<BindingScan>(); }
        return loadScans(g, genomeMap.getScans(g.getDBID()));
    }
    
    public Collection<BindingScan> loadScans(Genome g, Collection<Integer> dbids) throws SQLException { 
        LinkedList<BindingScan> scans = new LinkedList<BindingScan>();
        for(int dbid : dbids) { scans.addLast(loadScan(g, dbid)); }
        return scans;
    }
    
    public Collection<BindingScan> loadScans(Genome g, String version, String type) throws SQLException { 
        LinkedList<BindingScan> scans = new LinkedList<BindingScan>();
        synchronized (loadScansByVersionType) {
            loadScansByVersionType.setString(1, version);
            loadScansByVersionType.setString(2, type);
            ResultSet rs = loadScansByVersionType.executeQuery();
            while(rs.next()) { 
                BindingScan bs = new BindingScan(g, rs);
                scans.addLast(bs);
            }
            rs.close();
        }
        return scans;
    }
    
    public Collection<BindingScan> searchScans(Genome g, String version, String type) throws SQLException { 
        LinkedList<BindingScan> scans = new LinkedList<BindingScan>();        
        PreparedStatement ps = BindingScan.prepareLoadByLikeVersionType(c);
        ps.setString(1, version);
        ps.setString(2, type);
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            BindingScan bs = new BindingScan(g, rs);
            scans.addLast(bs);
        }
        rs.close();
        ps.close();
        return scans;
    }
    
	public BindingScan loadScan(Genome genome, int dbid) throws SQLException {
        //if(!genomeMap.containsScan(dbid)) { throw new IllegalArgumentException(); }
        
        // This is commented out, because we don't *really* want to check that there are *already*
        // regions for the given genome+scan in the database: when we want to scan a *new* genome
        // for a particular scan, we'll have to load the BindingScan in that new genome first, and it 
        // won't be in the genomeMap until we actually enter a scan for it in the new genome.
        // Think about it.
        //if(!genomeMap.getGenomes(dbid).contains(genome.getDBID())) { throw new IllegalArgumentException(); }
        BindingScan bs;
        synchronized (loadScanByID) {
            loadScanByID.setInt(1, dbid);
            ResultSet rs = loadScanByID.executeQuery();
            if(!rs.next()) { throw new IllegalArgumentException(String.valueOf(dbid)); }
            bs = new BindingScan(genome, rs);
            rs.close();
        }
		return bs;
	}
    
	public Collection<BindingScan> loadScans(Region r) throws SQLException {
        Genome genome = r.getGenome();
        
        LinkedList<BindingScan> scans = new LinkedList<BindingScan>();
        synchronized (loadScanByRegion) {
            loadScanByRegion.setInt(1, genome.getChromID(r.getChrom()));
            loadScanByRegion.setInt(2, r.getStart());
            loadScanByRegion.setInt(3, r.getStart());
            loadScanByRegion.setInt(4, r.getStart());
            loadScanByRegion.setInt(5, r.getEnd());        	   
            ResultSet rs = loadScanByRegion.executeQuery();
            while(rs.next()) { 
                BindingScan bs = new BindingScan(genome, rs);
                scans.addLast(bs);
            }
            rs.close();
        }
        
		return scans;
	}
    
    public Collection<Pair<Integer,Integer>> loadTypedExptPairs(BindingScan bs) throws SQLException { 
        Vector<Pair<Integer,Integer>> locs = new Vector<Pair<Integer,Integer>>();
        synchronized (loadExpts) {
            loadExpts.setInt(1, bs.getDBID());
            ResultSet rs = loadExpts.executeQuery();
            while(rs.next()) { 
                int type = rs.getInt(1);
                int expt = rs.getInt(2);
                locs.add(new Pair<Integer,Integer>(type, expt));
            }
            rs.close();
        }
        return locs;        
    }
    
    public Collection<ExptLocator> loadExpts(BindingScan bs) throws SQLException { 
        HashSet<ExptLocator> locs = new HashSet<ExptLocator>();
        synchronized (loadExpts) {
            loadExpts.setInt(1, bs.getDBID());
            ResultSet rs = loadExpts.executeQuery();
            while(rs.next()) { 
                int type = rs.getInt(1);
                int expt = rs.getInt(2);
                ExptLocator loc = lookupLocatorID(bs.getGenome(), type, expt);
                locs.add(loc);
            }
            rs.close();
        }
        return locs;
    }
    
    public Map<String,String> loadParams(BindingScan bs) throws SQLException {
        if(bs.getDBID() == -1) { throw new IllegalArgumentException(); }
        Map<String,String> params = new HashMap<String,String>();

        synchronized (loadParams) {
            loadParams.setInt(1, bs.getDBID());
            ResultSet rs = loadParams.executeQuery();
            while(rs.next()) { 
                String key = rs.getString(1);
                String value = rs.getString(2);
                params.put(key, value);
            }
            rs.close();
        }
        
        return params;
    }
    
    public Collection<Region> loadRegions(BindingScan bs) throws SQLException { 
        if(bs.getDBID() == -1) { throw new IllegalArgumentException(); }
        LinkedList<Region> regions = new LinkedList<Region>();
        
        synchronized(loadRegions) {
            loadRegions.setInt(1, bs.getDBID());
            ResultSet rs = loadRegions.executeQuery();
            while(rs.next()) { 
                int chromID = rs.getInt(1);
                int start = rs.getInt(2);
                int end = rs.getInt(3);
                String chromName = bs.getGenome().getChromName(chromID);
                Region r = new Region(bs.getGenome(), chromName, start, end);
                regions.addLast(r);
            }
            rs.close();
        }
        
        return regions;
    }
	
	public Collection<BindingEvent> loadEvents(BindingScan bs) throws SQLException {
        LinkedList<BindingEvent> events = new LinkedList<BindingEvent>();

        Genome genome = bs.getGenome();
		String type = bs.getVersion();
        synchronized (loadEventByScan) {
            loadEventByScan.setInt(1, bs.getDBID());
            ResultSet rs = loadEventByScan.executeQuery();
            while(rs.next()) { 
                int chromID = rs.getInt(2);
                String chrom = genome.getChromName(chromID);
                int start = rs.getInt(3);
                int end = rs.getInt(4);
                double size = rs.getDouble(5);
                double conf = rs.getDouble(6);
                BindingEvent evt = new BindingEvent(genome, chrom, start, end, size, conf, type);
                events.addLast(evt);
            }
            rs.close();
        }
		return events;
	}
    
    public Collection<BindingEvent> loadEvents(BindingScan bs, Region r) throws SQLException {
        if(!r.getGenome().equals(bs.getGenome())) { throw new IllegalArgumentException(); }
        
        LinkedList<BindingEvent> events = new LinkedList<BindingEvent>();
        
        Genome genome = bs.getGenome();
        synchronized (loadEventByScanRegion) {
            loadEventByScanRegion.setInt(1, bs.getDBID());
            loadEventByScanRegion.setInt(2, genome.getChromID(r.getChrom()));
            loadEventByScanRegion.setInt(3, r.getStart());
            loadEventByScanRegion.setInt(4, r.getStart());
            loadEventByScanRegion.setInt(5, r.getStart());
            loadEventByScanRegion.setInt(6, r.getEnd());
        
            String type = bs.getVersion();
            ResultSet rs = loadEventByScanRegion.executeQuery();
            while(rs.next()) { 
                int chromID = rs.getInt(2);
                String chrom = genome.getChromName(chromID);
                int start = rs.getInt(3);
                int end = rs.getInt(4);
                double size = rs.getDouble(5);
                double conf = rs.getDouble(6);
                BindingEvent evt = new BindingEvent(genome, chrom, start, end, size, conf, type);
                events.addLast(evt);
            }
            rs.close();
        }

        return events;
    }
	
	public void insertScan(BindingScan bs) throws SQLException { 
		bs.insertIntoDB(c, genomeMap, insertScan);
	}
    
    public void insertNewRegion(BindingScan bs, Region r) throws SQLException { 
        int dbid = bs.getDBID();
        if(dbid == -1) { throw new IllegalArgumentException(); }
        synchronized (insertRegion) {
            insertRegion.setInt(1, dbid);
            insertRegion.setInt(2, r.getGenome().getChromID(r.getChrom()));
            insertRegion.setInt(3, r.getStart());
            insertRegion.setInt(4, r.getEnd());
            
            insertRegion.executeUpdate();
        }
        
        genomeMap.insertNewMapping(c, bs.getDBID(), r.getGenome().getDBID());
    }
    
    public void insertNewRegions(BindingScan bs, Collection<Region> regions) throws SQLException { 
        for(Region r : regions) { insertNewRegion(bs, r); }
    }
    
    public void insertNewParam(BindingScan bs, String k, String v) throws SQLException { 
        int dbid = bs.getDBID();
        if(dbid == -1) { throw new IllegalArgumentException(); }
    
        synchronized (insertParam) {
            insertParam.setInt(1, dbid);
            insertParam.setString(2, k);
            insertParam.setString(3, v.length() > 49 ? v.substring(v.length()-49, v.length()) : v);
            
            insertParam.executeUpdate();
        }
    }
    
    public void insertNewExpts(BindingScan bs, Collection<ExptLocator> locs) throws SQLException { 
        for(ExptLocator loc : locs) { insertNewExpt(bs, loc); }
    }
    
    public void insertNewExpt(BindingScan bs, ExptLocator loc) throws SQLException { 
        int type = BindingScan.getLocatorType(loc);
        int[] expts = BindingScan.getLocatorIDs(loc, c);
        for(int i = 0; i < expts.length; i++) { 
            insertNewExpt(bs, expts[i], type);
        }
    }
    
    public void insertNewExpts(BindingScan bs, int[] expts, int type) throws SQLException { 
        for(int i = 0; i < expts.length; i++) { 
            insertNewExpt(bs, expts[i], type); 
        }
    }
    
    public void insertNewExpt(BindingScan bs, int expt, int type) throws SQLException { 
        int dbid = bs.getDBID();
        if(dbid == -1) { throw new IllegalArgumentException(); }

        synchronized (insertExpt) {
            insertExpt.setInt(1, dbid);
            insertExpt.setInt(2, expt);
            insertExpt.setInt(3, type);
            insertExpt.executeUpdate();
        }
    }
    
    public void insertNewParams(BindingScan bs, Map<String,String> params) throws SQLException { 
        for(String k : params.keySet()) { insertNewParam(bs, k, params.get(k)); }
    }
    
    public void insertExpt(BindingScan bs, ExptLocator loc) throws SQLException { 
        int type = BindingScan.getLocatorType(loc);
        int[] expts = BindingScan.getLocatorIDs(loc, c);
    }
    
    public void insertEvent(BindingScan bs, BindingEvent evt) throws SQLException { 
        if(bs.getDBID() == -1) { throw new IllegalArgumentException(); }
        if(!evt.getGenome().equals(bs.getGenome())) { throw new IllegalArgumentException(); }
        
        synchronized (insertEvent) {
            insertEvent.setInt(1, bs.getDBID());
            insertEvent.setInt(2, evt.getGenome().getChromID(evt.getChrom()));
            insertEvent.setInt(3, evt.getStart());
            insertEvent.setInt(4, evt.getEnd());
            insertEvent.setDouble(5, evt.getSize());
            insertEvent.setDouble(6, evt.getConf());
            
            insertEvent.executeUpdate();
        }
    }
    
    public void deleteScan(int dbid) throws SQLException {
        c.setAutoCommit(false);
        
        Statement s = c.createStatement();
        
        s.executeUpdate("delete from bindingscanparam where scan=" + dbid);
        s.executeUpdate("delete from bindingscanToExpt where scan=" + dbid);
        s.executeUpdate("delete from bindingscanregion where scan=" + dbid);
        s.executeUpdate("delete from bindingevent where scan=" + dbid);
        s.executeUpdate("delete from bindingscanToGenome where scan=" + dbid);
        s.executeUpdate("delete from bindingscan where id=" + dbid);
        
        s.close();
        c.commit();
        genomeMap.removeScan(dbid);
        
        c.setAutoCommit(true);
    }
    
    public ExptLocator lookupLocatorID(Genome g, int type, int id) throws SQLException {
        ExptLocator loc = null;
        
        Statement s = c.createStatement();
        ResultSet rs;
        String name, version;

        switch(type) { 
        case BindingScan.AGILENT_TYPE:
            rs = s.executeQuery("select name, version from experiment where id=" + id);
            if(rs.next()) { 
                name = rs.getString(1);
                version = rs.getString(2);
                loc = new ChipChipLocator(g, name, version);
            }
            rs.close();
            break;
        case BindingScan.BAYES_TYPE:
            rs = s.executeQuery("select name, version from bayesanalysis where id=" + id);
            if(rs.next()) { 
                name = rs.getString(1);
                version = rs.getString(2);
                loc = new BayesLocator(g, name, version);
            }
            rs.close();
            break;
        case BindingScan.MSP_TYPE:
            rs = s.executeQuery("select name, version from rosettaanalysis where id=" + id);
            if(rs.next()) { 
                name = rs.getString(1);
                version = rs.getString(2);
                loc = new MSPLocator(g, name, version);
            }
            rs.close();
            break;
        }
        
        s.close();
        
        return loc;
    }
	
	public static PreparedStatement prepareInsertEvent(java.sql.Connection c) throws SQLException { 
		String str = "insert into bindingevent (scan, chromosome, startpos, stoppos, eventsize, eventconf) " + 
			" values (?, ?, ?, ?, ?, ?)";
		return c.prepareStatement(str);
	}
	
	public static PreparedStatement prepareLoadEventByScan(java.sql.Connection c) 
	throws SQLException { 
	    String query = "select scan, chromosome, startpos, stoppos, eventsize, eventconf from bindingevent " +
	    "where scan=? order by chromosome,startpos";
	    return c.prepareStatement(query);
	}
	
	public static PreparedStatement prepareLoadEventByScanRegion(java.sql.Connection c) 
	throws SQLException { 
	    String query = "select scan, chromosome, startpos, stoppos, eventsize, eventconf from bindingevent " +
	    "where scan=? and chromosome=? and ((startpos <= ? and stoppos >= ?) or (startpos >=? and startpos <=?)) " +
        "order by chromosome,startpos";
	    return c.prepareStatement(query);
	}
    
    public static PreparedStatement loadExptsByScan(java.sql.Connection c) throws SQLException {
        String query = "select scanexpt, scantype from bindingscanToExpt where scan=?";
        return c.prepareStatement(query);
    }
    
    public static PreparedStatement loadParamsByScan(java.sql.Connection c) throws SQLException {
        String query = "select key, value from bindingscanparam where scan=?";
        return c.prepareStatement(query);
    }
    
    public static PreparedStatement loadRegionsByScan(java.sql.Connection c) throws SQLException { 
        String query = "select chromosome, startpos, stoppos from bindingscanregion where scan=? order by chromosome,startpos";
        return c.prepareStatement(query);
    }
}
