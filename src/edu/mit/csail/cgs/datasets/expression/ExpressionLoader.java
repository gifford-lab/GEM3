/*
 * Created on Mar 14, 2007
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.TimeSeriesLoader;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class ExpressionLoader implements Closeable {
    
    public static void main(String[] args) { 
        try {
            ExpressionLoader loader = new ExpressionLoader();
            
            System.out.println("Platforms:");
            Collection<ProbePlatform> platforms = loader.loadProbePlatforms();
            for(ProbePlatform platform : platforms) { 
                System.out.println("\t\"" + platform.getName() + "\" (" + platform.getType() + ")");
            }
            System.out.println();
            
            System.out.println("Experiments:");
            Collection<Experiment> expts = loader.loadAllExperiments();
            for(Experiment expt : expts) { 
                System.out.println("\t\"" + expt.getName() + "\" --> \"" + expt.getPlatform().getName() + "\"");
            }
            System.out.println();
            
            loader.close();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public static final String role = "expression";
    
    private MetadataLoader chipLoader;
    private TimeSeriesLoader timeLoader;

    private boolean ownsConnection;
    private java.sql.Connection cxn;
    
    private PreparedStatement loadExptByID; 
    private PreparedStatement loadExptParamsByID; 
    private PreparedStatement loadProbeByID, loadProbesByExptID; 
    private PreparedStatement loadMeasureByBothIDs, loadMeasureByExptID, loadMeasureByProbeID; 
    private PreparedStatement loadProcessingByID; 
    private PreparedStatement loadProcessingParamsByID; 
    private PreparedStatement loadProcessingInputsByID; 
    private PreparedStatement loadSetByID, loadSetMembersByID;
    private PreparedStatement loadPlatformByID;
    
    private Map<Integer,Experiment> cachedExpts;
    private Map<Integer, Processing> cachedProcessings;
    private Map<Integer,ExperimentSet> cachedSets;
    private Map<Integer,ProbePlatform> cachedPlatforms;
    
    public ExpressionLoader(java.sql.Connection cxn) throws SQLException { 
        chipLoader = new MetadataLoader();
        timeLoader = new TimeSeriesLoader();
        this.cxn = cxn;
        ownsConnection = false;
        init();
    }
    
    public ExpressionLoader() throws SQLException { 
        chipLoader = new MetadataLoader();
        timeLoader = new TimeSeriesLoader();
        
        try {
            cxn = DatabaseFactory.getConnection(role);
            ownsConnection = true;
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
        
        init();
    }
    
    public java.sql.Connection getConnection() { return cxn; }
    
    private void init() throws SQLException { 
        loadExptByID = Experiment.prepareLoadByID(cxn);
        loadExptParamsByID = Experiment.prepareLoadParamsByID(cxn);
        
        loadProbeByID = Probe.prepareLoadByID(cxn);
        loadProbesByExptID = Probe.prepareLoadByExperimentID(cxn);
        
        loadMeasureByBothIDs = ExprMeasurement.prepareLoadByBothIDs(cxn);
        loadMeasureByExptID = ExprMeasurement.prepareLoadByExptID(cxn);
        loadMeasureByProbeID = ExprMeasurement.prepareLoadByProbeID(cxn);
        
        loadProcessingByID = Processing.prepareLoadByID(cxn);
        loadProcessingParamsByID = Processing.prepareLoadParamsByID(cxn);
        loadProcessingInputsByID = Processing.prepareLoadInputsByID(cxn);
        
        loadSetByID = ExperimentSet.prepareLoadByID(cxn);
        loadSetMembersByID = ExperimentSet.prepareLoadMembersByID(cxn);
        
        loadPlatformByID = ProbePlatform.prepareLoadByID(cxn);
        
        cachedExpts = new HashMap<Integer,Experiment>();
        cachedProcessings = new HashMap<Integer, Processing>();
        cachedSets = new HashMap<Integer,ExperimentSet>();
        cachedPlatforms = new HashMap<Integer,ProbePlatform>();        
    }
    
    public void close() {
    	
    	cachedExpts.clear();
    	cachedProcessings.clear();
    	cachedSets.clear();
    	cachedPlatforms.clear();
        
        chipLoader.close();
        chipLoader = null;
        
        timeLoader.close();
        timeLoader = null;
        
        try {
            loadExptByID.close(); loadExptByID = null;
            loadExptParamsByID.close();loadExptParamsByID = null;
        
            loadProbeByID.close(); loadProbeByID = null;
            loadProbesByExptID.close(); loadProbesByExptID = null;
    
            loadMeasureByBothIDs.close(); loadMeasureByBothIDs = null;
            loadMeasureByExptID.close(); loadMeasureByExptID = null;
            loadMeasureByProbeID.close(); loadMeasureByProbeID = null;
            
            loadProcessingByID.close(); loadProcessingByID = null;
            loadProcessingParamsByID.close(); loadProcessingParamsByID = null;
            loadProcessingInputsByID.close(); loadProcessingInputsByID = null;
            
            loadSetByID.close(); loadSetByID = null;
            loadSetMembersByID.close(); loadSetMembersByID = null;
            
            loadPlatformByID.close(); loadPlatformByID = null;
            
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
        if(ownsConnection) { 
            DatabaseFactory.freeConnection(cxn);
        }
        
        cxn = null;
    }
    
    public boolean isClosed() {
        return cxn == null;
    }
    
    public Processing loadProcessing(int dbid) throws SQLException {
    	if (cachedProcessings.containsKey(dbid)) {return cachedProcessings.get(dbid); }
    	
    	Processing p = null;

        synchronized (loadProcessingByID) {
            loadProcessingByID.setInt(1, dbid);
            loadProcessingParamsByID.setInt(1, dbid);
            loadProcessingInputsByID.setInt(1, dbid);
            
            ResultSet rs = loadProcessingByID.executeQuery();
            ResultSet prs = loadProcessingParamsByID.executeQuery();
            ResultSet irs = loadProcessingInputsByID.executeQuery();
            
            if (!rs.next()) {
                throw new IllegalArgumentException("Processing ID: " + dbid);
            } else {
                p = new Processing(rs, prs, irs, this); 
                cachedProcessings.put(dbid, p);
            }
            
            rs.close();
            prs.close();
            irs.close();
        }
    	
    	return p;
    }
    
    public ProbePlatform loadPlatform(int dbid) throws SQLException {
    	if(cachedPlatforms.containsKey(dbid)) { return cachedPlatforms.get(dbid); }
    	ProbePlatform plat = null;

        synchronized (loadPlatformByID) {
            loadPlatformByID.setInt(1, dbid);
            ResultSet rs = loadPlatformByID.executeQuery();
            
            if(rs.next()) { 
                plat = new ProbePlatform(rs);
            }
            
            rs.close();
        }
    	if(plat == null) { throw new IllegalArgumentException("Unknwon Platform : " + dbid); }
    	return plat;
    }
    
    public ProbePlatform loadPlatform(String name) throws SQLException {
        ProbePlatform platform = null;
        PreparedStatement ps = cxn.prepareStatement("select id, name, type from probe_platform where name=?");
        ps.setString(1, name);
        ResultSet rs = ps.executeQuery();
        
        if(rs.next()) { 
            int id = rs.getInt(1);
            if(cachedPlatforms.containsKey(id)) { 
                platform = cachedPlatforms.get(id);
            } else { 
                platform = new ProbePlatform(rs);
                cachedPlatforms.put(id, platform);
            }
        }
        
        rs.close();
        ps.close();
        
        return platform;
    }
    
    public ExperimentSet loadExperimentSet(int dbid) throws SQLException { 
        if(cachedSets.containsKey(dbid)) { return cachedSets.get(dbid); }
        
        ExperimentSet set = null;
        synchronized (loadSetByID) {
            loadSetByID.setInt(1, dbid);
            loadSetMembersByID.setInt(1, dbid);
            
            ResultSet rs = loadSetByID.executeQuery();
            ResultSet mrs = loadSetMembersByID.executeQuery();
            
            if(rs.next()) { 
                set = new ExperimentSet(rs, mrs, this);
            }
            
            rs.close();
            mrs.close();
        }
        
        if(set == null) { throw new IllegalArgumentException("Unknown ExperimentSet dbid: " + dbid); }
        return set;
    }
    
    public Map<String,Object> loadExperimentStatistics(Experiment expt) throws SQLException { 
        Map<String,Object> stats = new HashMap<String,Object>();
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select max(value), min(value) from " +
                "measurement where experiment=" + expt.getDBID());
        rs.next();
        double max = rs.getDouble(1);
        double min = rs.getDouble(2);
        rs.close();
        
        stats.put("max", max);
        stats.put("min", min);
        
        s.close();
        return stats;
    }

	public Experiment loadExperiment(String name) throws SQLException { 
		Experiment expt = null;
    	String query = "select id from experiment where name=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, name);
		ResultSet rs = ps.executeQuery();

		if(rs.next()) { 
			int id = rs.getInt(1);
			expt = loadExperiment(id);
		}

		rs.close();
		ps.close();
		return expt;
	}
    
    public Experiment loadExperiment(int dbid) throws SQLException { 
        if(cachedExpts.containsKey(dbid)) { return cachedExpts.get(dbid); }
        
        Experiment e = null;
        synchronized (loadExptByID) {
            loadExptByID.setInt(1, dbid);
            loadExptParamsByID.setInt(1, dbid);
            
            ResultSet rs = loadExptByID.executeQuery();
            ResultSet prs = loadExptParamsByID.executeQuery();
            
            if(!rs.next()) { 
                System.err.println("Couldn't find experiment with ID: " + dbid);
            } else { 
                e = new Experiment(rs, prs, this, timeLoader, chipLoader);
                cachedExpts.put(dbid, e);
            }
            
            rs.close();
            prs.close();
        }
        
        return e;
    }
    
    public Collection<Experiment> loadExperiments(Collection<Integer> dbids) throws SQLException { 
    	LinkedList<Experiment> expts = new LinkedList<Experiment>();
    	for(int dbid : dbids) { 
    		Experiment e = loadExperiment(dbid);
    		if(e != null) { expts.addLast(e); }
    	}
    	return expts;
    }
    
    public Probe loadProbe(int dbid) throws SQLException { 
        Probe p = null;
        synchronized (loadProbeByID) {
            loadProbeByID.setInt(1, dbid);
            ResultSet rs = loadProbeByID.executeQuery();
            if(!rs.next()) { 
                throw new IllegalArgumentException("Unknown Probe DBID: " + dbid);
            } else { 
                p = new Probe(rs, this);
            }
            rs.close();
        }
        return p;
    }
    
    public Collection<LocatedExprMeasurement> loadMeasurementsInRegion(Region r, ProbePlatform pp, Experiment expt) 
        throws SQLException { 

        LinkedList<LocatedExprMeasurement> measurements = new LinkedList<LocatedExprMeasurement>();
        PreparedStatement ps = cxn.prepareStatement("select p.id, p.name, p.platform, " +
                "pl.startpos, pl.stoppos, em.value from " +
                "probe p, probe_location pl, measurement em where em.probe=p.id and p.id=pl.probe " +
                "and em.experiment=? and p.platform=? " +
                "and pl.chromosome=? and " +
                "((pl.startpos <= ? and pl.stoppos >= ?) or " +
                "(pl.startpos >= ? and pl.startpos <= ?))");
        
        ps.setInt(1, expt.getDBID());
        ps.setInt(2, pp.getDBID());
        ps.setInt(3, r.getGenome().getChromID(r.getChrom()));
        ps.setInt(4, r.getStart());
        ps.setInt(5, r.getStart());
        ps.setInt(6, r.getStart());
        ps.setInt(7, r.getEnd());
        
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            Probe probe = new Probe(rs, this);
            int start = rs.getInt(4), end = rs.getInt(5);
            double value = rs.getDouble(6);
            NamedRegion region = 
                new NamedRegion(r.getGenome(), r.getChrom(), start, end, probe.getName());
            LocatedExprMeasurement lem = new LocatedExprMeasurement(expt, probe, value, region);
            measurements.addLast(lem);
        }
        rs.close();
                
        ps.close();
        return measurements;
    }
    
    public Collection<LocatedProbe> loadProbesInRegion(Region r, ProbePlatform pp) throws SQLException {
    	LinkedList<LocatedProbe> probes = new LinkedList<LocatedProbe>();
    	PreparedStatement ps = cxn.prepareStatement("select p.id, p.name, p.platform, " +
    			"pl.startpos, pl.stoppos from " +
    			"probe p, probe_location pl where p.platform=? and p.id=pl.probe " +
    			"and pl.chromosome=? and " +
    			"((pl.startpos <= ? and pl.stoppos >= ?) or " +
    			"(pl.startpos >= ? and pl.startpos <= ?))");
    	
    	ps.setInt(1, pp.getDBID());
    	ps.setInt(2, r.getGenome().getChromID(r.getChrom()));
    	ps.setInt(3, r.getStart());
    	ps.setInt(4, r.getStart());
    	ps.setInt(5, r.getStart());
    	ps.setInt(6, r.getEnd());
    	
    	ResultSet rs = ps.executeQuery();
    	while(rs.next()) { 
    		Probe probe = new Probe(rs, this);
    		int start = rs.getInt(4), end = rs.getInt(5);
    		NamedRegion region = 
    			new NamedRegion(r.getGenome(), r.getChrom(), start, end, probe.getName());
    		probes.addLast(new LocatedProbe(probe, region));
    	}
    	rs.close();
    			
    	ps.close();
    	return probes;
    }
    
    public Collection<NamedRegion> loadProbeLocations(Probe p, Genome g) throws SQLException { 
    	LinkedList<NamedRegion> regions = new LinkedList<NamedRegion>();
    	PreparedStatement ps = cxn.prepareStatement("select chromosome, startpos, stoppos from probe_location where probe=?");

    	ps.setInt(1, p.getDBID());
    	ResultSet rs = ps.executeQuery();
    	while(rs.next()) { 
    		int chromID = rs.getInt(1);
    		String chrom = g.getChromName(chromID);
    		int start = rs.getInt(2), end = rs.getInt(3);
    		NamedRegion r = new NamedRegion(g, chrom, start, end, p.getName());
    		regions.addLast(r);
    	}
    	rs.close();
    	
    	ps.close();
    	return regions;
    }
    
    public Map<Probe,Set<NamedRegion>> loadAllProbeLocations(Set<Probe> pps, Genome g) throws SQLException { 
    	Map<Probe,Set<NamedRegion>> map = new HashMap<Probe,Set<NamedRegion>>();
    	PreparedStatement ps = cxn.prepareStatement("select chromosome, startpos, stoppos from probe_location where probe=?");

    	for(Probe p : pps) { 
        	Set<NamedRegion> regions = new HashSet<NamedRegion>();
    		ps.setInt(1, p.getDBID());
    		ResultSet rs = ps.executeQuery();
    		while(rs.next()) { 
    			int chromID = rs.getInt(1);
    			String chrom = g.getChromName(chromID);
    			int start = rs.getInt(2), end = rs.getInt(3);
    			NamedRegion r = new NamedRegion(g, chrom, start, end, p.getName());
    			regions.add(r);
    		}
    		rs.close();
    		
    		map.put(p, regions);
    	}
    	
    	ps.close();
    	return map;
    }
    
    public Collection<Experiment> loadAllExperiments(Probe p) throws SQLException { 
    	LinkedList<Experiment> exptList = new LinkedList<Experiment>();
    	String query = "select e.id from experiment e, probe_platform pp, " +
    			"probe p where p.id=? and p.platform=pp.id and e.platform=pp.id";
    	PreparedStatement ps = cxn.prepareStatement(query);
    	ps.setInt(1, p.getDBID());
    	ResultSet rs = ps.executeQuery();
    	while(rs.next()) { 
    		int id = rs.getInt(1);
    		Experiment expt = loadExperiment(id);
    		exptList.addLast(expt);
    	}
    	rs.close();
    	
    	ps.close();
    	return exptList;
    }
    
    public Collection<Experiment> loadAllExperiments() throws SQLException {
        LinkedList<Experiment> exptList = new LinkedList<Experiment>();
        String query = "select id from experiment order by name";
        PreparedStatement ps = cxn.prepareStatement(query);
        
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            int id = rs.getInt(1);
            exptList.addLast(loadExperiment(id));
        }
        rs.close();
        
        ps.close();
        return exptList;
    }
    
    public void deleteExperiment(int dbid) throws SQLException {
        cxn.setAutoCommit(false);
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select distinct(ip.processing) from processing_inputpair ip " +
                "where (ip.input=" + dbid + " or ip.output=" + dbid + ")");
        Set<Integer> procids = new HashSet<Integer>();
        while(rs.next()) { procids.add(rs.getInt(1)); }
        rs.close();

        s.executeUpdate("delete from processing_inputpair where (input=" + dbid + " or output=" + dbid + ")");
        for(int proc : procids) { 
            s.executeUpdate("delete from processing_params where processing=" + proc);
            s.executeUpdate("delete from processing where id=" + proc);
        }
        
        s.executeUpdate("delete from measurement where experiment=" + dbid);
        s.executeUpdate("delete from experiment_params where experiment=" + dbid);
        s.executeUpdate("delete from experiment_set_member where experiment=" + dbid);
        s.executeUpdate("delete from experiment where id=" + dbid);
        
        if(cachedExpts.containsKey(dbid)) { cachedExpts.remove(dbid); }
        
        s.close();
        cxn.commit();
        cxn.setAutoCommit(true);
    }
    
    public void deleteExperimentSet(int dbid) throws SQLException { 
        cxn.setAutoCommit(false);
        Statement s = cxn.createStatement();
        
        s.executeUpdate("delete from experiment_set_member where experiment_set=" + dbid);
        s.executeUpdate("delete from experiment_set where id=" + dbid);
        
        if(cachedSets.containsKey(dbid)) { cachedSets.remove(dbid); } 
        
        s.close();
        cxn.commit();
        cxn.setAutoCommit(true);
    }
    
    public void deleteProbePlatform(int dbid) throws SQLException { 
        cxn.setAutoCommit(false);
        Statement s = cxn.createStatement();

        s.executeUpdate("delete from probe where platform=" + dbid);
        s.executeUpdate("delete from probe_platform where id=" + dbid);
        
        if(cachedPlatforms.containsKey(dbid)) { cachedPlatforms.remove(dbid); } 
        
        s.close();
        cxn.commit();
        cxn.setAutoCommit(true);        
    }
    
    public Collection<Probe> loadProbesByName(String name) throws SQLException { 
    	LinkedList<Probe> pbs = new LinkedList<Probe>();
        PreparedStatement s = cxn.prepareStatement("select id, name, platform from probe where name=?");
        s.setString(1, name);
        ResultSet rs = s.executeQuery();
        if(rs.next()) { 
        	Probe p = new Probe(rs, this);
        	pbs.addLast(p);
        }
        rs.close();
        s.close();    	
    	return pbs;
    }
    
    public Map<String,Probe> loadProbes(Set<String> probeNames) throws SQLException { 
        HashMap<String,Probe> probemap = new HashMap<String,Probe>();
        
        PreparedStatement s = cxn.prepareStatement("select id, name, platform from probe where name=?");
        for(String name : probeNames) { 
            s.setString(1, name);
            ResultSet rs = s.executeQuery();
            if(rs.next()) { 
                Probe p = new Probe(rs, this);
                probemap.put(name, p);
            }
            rs.close();
        }
        s.close();
        
        return probemap;
    }
    
    public Map<String,Probe> loadProbes(ProbePlatform platform, Set<String> probeNames) throws SQLException { 
        HashMap<String,Probe> probemap = new HashMap<String,Probe>();
        
        PreparedStatement s = cxn.prepareStatement("select id, name, platform from probe where platform=? and name=?");
        for(String name : probeNames) { 
            s.setInt(1, platform.getDBID());
            s.setString(2, name);

            ResultSet rs = s.executeQuery();
            if(rs.next()) { 
                Probe p = new Probe(rs, this);
                probemap.put(name, p);
            }
            rs.close();
        }
        s.close();
        
        return probemap;
    }
    
    public Collection<ProbePlatform> loadProbePlatforms() throws SQLException {
        LinkedList<ProbePlatform> plats = new LinkedList<ProbePlatform>();
        PreparedStatement ps = cxn.prepareStatement("select id, name, type from probe_platform order by name");
        ResultSet rs = ps.executeQuery();
        
        while(rs.next()) { 
            int id = rs.getInt(1);
            if(cachedPlatforms.containsKey(id)) { 
                plats.addLast(cachedPlatforms.get(id));
            } else { 
                ProbePlatform pp = new ProbePlatform(rs);
                cachedPlatforms.put(id, pp);
                plats.addLast(pp);
            }
        }
        
        rs.close();
        ps.close();
        return plats;
    }
    
    public Collection<Processing> loadProcessingWithInputExperiment(Experiment e) throws SQLException { 
    	LinkedList<Processing> procs = new LinkedList<Processing>();
    	
    	Statement s = cxn.createStatement();
    	ResultSet rs = s.executeQuery("select distinct(proc.id) from processing proc, " +
    			"processing_inputpair ip where proc.id=ip.processing and ip.input=" + e.getDBID());
    	
    	while(rs.next()) { 
    		int dbid = rs.getInt(1);
    		Processing proc = loadProcessing(dbid);
    		procs.addLast(proc);
    	}
    	
    	rs.close();
    	s.close();
    	
    	return procs;
    }

    public Collection<Processing> loadProcessingWithOutputExperiment(Experiment e) throws SQLException { 
    	LinkedList<Processing> procs = new LinkedList<Processing>();
    	
    	Statement s = cxn.createStatement();
    	ResultSet rs = s.executeQuery("select distinct(proc.id) from processing proc, " +
    			"processing_inputpair ip where proc.id=ip.processing and ip.output=" + e.getDBID());
    	
    	while(rs.next()) { 
    		int dbid = rs.getInt(1);
    		Processing proc = loadProcessing(dbid);
    		procs.addLast(proc);
    	}
    	
    	rs.close();
    	s.close();
    	
    	return procs;
    }

    /**
     * Returns a collection of Experiment objects in which this Probe was measured.
     * 
     * @param p
     * @return
     * @throws SQLException
     */
    public Collection<Experiment> loadProbeExperiments(Probe p) throws SQLException {
        LinkedList<Experiment> expts = new LinkedList<Experiment>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select m.experiment from measurement m where m.probe=" + p.getDBID());
        while(rs.next()) { 
            int exptID = rs.getInt(1);
            Experiment e = loadExperiment(exptID);
            expts.addLast(e);
        }
        rs.close();
        s.close();
        return expts;
    }
    
    /**
     * Returns an ExprMeasurement object corresponding to the given Experiment and Probe objects.
     * 
     * @param expt
     * @param p
     * @return
     * @throws SQLException
     */
    public ExprMeasurement loadMeasurement(Experiment expt, Probe p) throws SQLException { 
        ExprMeasurement em = null;
        synchronized (loadMeasureByBothIDs) {
            loadMeasureByBothIDs.setInt(1, expt.getDBID());
            loadMeasureByBothIDs.setInt(2, p.getDBID());
            
            ResultSet rs = loadMeasureByBothIDs.executeQuery();
            if(rs.next()) { 
                em = new ExprMeasurement(expt, p, rs.getDouble(3));
            }
            rs.close();
        }
        return em;
    }

    /**
     * Loads an ExprMeasurement object for *each* measurement made in a particular Experiment.
     * 
     * @param expt
     * @return
     * @throws SQLException
     */
    public Collection<ExprMeasurement> loadExperimentMeasurements(Experiment expt) throws SQLException { 

        LinkedList<ExprMeasurement> ms = new LinkedList<ExprMeasurement>();

        synchronized (loadMeasureByExptID) {
            loadMeasureByExptID.setInt(1, expt.getDBID());
            ResultSet rs = loadMeasureByExptID.executeQuery();
            
            while(rs.next()) {
                // This is a little tricky.  See the comment above ExprMeasurement.prepareLoadByExptID(),
                // for a bit more explanation.  I'm doing it this way to avoid having to make *two* database 
                // hits for each particular measurement we want to read.
                
                Probe p = new Probe(rs, this);
                double value = rs.getDouble(4);
                
                ExprMeasurement em = new ExprMeasurement(expt, p, value);
                ms.addLast(em);
            }
            rs.close();
        }
        return ms;
    }
    
    /**
     * Loads *all* measured Probe objects, for the given Experiment.
     * 
     * @param expt
     * @return
     * @throws SQLException
     */
    public Collection<Probe> loadExperimentProbes(Experiment expt) throws SQLException { 
        LinkedList<Probe> probes = new LinkedList<Probe>();
        
        synchronized (loadProbesByExptID) {
            loadProbesByExptID.setInt(1, expt.getDBID());
            ResultSet rs = loadProbesByExptID.executeQuery();
            while(rs.next()) { 
                Probe p = new Probe(rs, this);
                probes.addLast(p);
            }
            rs.close();
        }
        return probes;
    }

}
