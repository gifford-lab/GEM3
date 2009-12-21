/*
 * Created on Mar 20, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.ArrayList;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class ChipChipMetadataLoader implements edu.mit.csail.cgs.utils.Closeable {

    public static final String role = "chipchip";

    private java.sql.Connection cxn;
    private MetadataLoader metaLoader;
    private boolean freeMetaLoader;
    
    private Map<Integer,ArrayDesign> cachedArrays;
    private Map<Integer,GALFile> cachedGALFiles;
    private Map<String,GALFile> nameCachedGALFiles;
    private Map<String,Set<ArrayDesign>> nameCachedArrays;

    public ChipChipMetadataLoader(MetadataLoader ml) throws SQLException {
        metaLoader = ml;
        freeMetaLoader = false;
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }
        
        cachedArrays = new HashMap<Integer,ArrayDesign>();
        cachedGALFiles = new HashMap<Integer,GALFile>();
        nameCachedGALFiles = new HashMap<String,GALFile>();
        nameCachedArrays = new HashMap<String,Set<ArrayDesign>>();
    }
    
    public ChipChipMetadataLoader() throws SQLException { 
        metaLoader = new MetadataLoader();
        freeMetaLoader = true;
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown role: " + role, e);
        }        

        cachedArrays = new HashMap<Integer,ArrayDesign>();
        cachedGALFiles = new HashMap<Integer,GALFile>();
        nameCachedGALFiles = new HashMap<String,GALFile>();
        nameCachedArrays = new HashMap<String,Set<ArrayDesign>>();
    }
    
    public void close() {
        if(freeMetaLoader) { metaLoader.close(); }
        metaLoader = null;
        
        cachedArrays.clear();
        cachedGALFiles.clear();
        nameCachedGALFiles.clear();
        nameCachedArrays.clear();
        
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }

    /* calls commit() on any underlying database connections to flush
       any database updates */
    public void commit() throws SQLException {
        cxn.commit();
    }

    public Collection<Cells> loadAllCells(Genome g) throws SQLException { 
    
        HashSet<Cells> values = new HashSet<Cells>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select e.cellsone, e.cellstwo from experiment e, " +
                                      "exptToGenome e2g where e.id=e2g.experiment and e2g.genome=" + g.getDBID()); 

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(metaLoader.loadCells(id));
            
            id = rs.getInt(2);
            values.add(metaLoader.loadCells(id));
        }

        rs.close();
        s.close();
        return values;
    }

    public Collection<Condition> loadAllConditions(Genome g) throws SQLException {
        
        HashSet<Condition> values = new HashSet<Condition>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select e.conditionone, e.conditiontwo from experiment e, " +
                                      "exptToGenome e2g where e.id=e2g.experiment and e2g.genome=" + g.getDBID()); 

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(metaLoader.loadCondition(id));
            
            id = rs.getInt(2);
            values.add(metaLoader.loadCondition(id));
        }

        rs.close();
        s.close();
        return values;
    }

    public Collection<Factor> loadAllFactors(Genome g) throws SQLException {
        
        HashSet<Factor> values = new HashSet<Factor>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select e.factorone, e.factortwo from experiment e, " +
                                      "exptToGenome e2g where e.id=e2g.experiment and e2g.genome=" + g.getDBID()); 

        while(rs.next()) { 
            int id = rs.getInt(1);
            values.add(metaLoader.loadFactor(id));

            id = rs.getInt(2);
            values.add(metaLoader.loadFactor(id));
        }

        rs.close();
        s.close();
        return values;
    }
    
    private void cacheArrayDesign(ArrayDesign ad) { 
        if(!cachedArrays.containsKey(ad.getDBID())) { 
            cachedArrays.put(ad.getDBID(), ad);
        }
        
        if(!nameCachedArrays.containsKey(ad.getName())) { 
            nameCachedArrays.put(ad.getName(), new HashSet<ArrayDesign>());
        }
        
        if(!nameCachedArrays.get(ad.getName()).contains(ad)) { 
            nameCachedArrays.get(ad.getName()).add(ad);
        }
    }
    
    private ArrayDesign retrieveArrayDesignFromCache(String name, Genome g) { 
        if(!nameCachedArrays.containsKey(name)) { return null; }
        for(ArrayDesign ad : nameCachedArrays.get(name)) { 
            if(ad.getGenome().equals(g)) { 
                return ad;
            }
        }
        return null;
    }

    public Collection<ArrayDesign> loadAllArrayDesigns(Organism spec) throws SQLException {

        HashSet<ArrayDesign> values = new HashSet<ArrayDesign>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id, name, genome from arraydesign");        
        Collection<Integer> ids = spec.getGenomeIDs();
        
        while (rs.next()) {
            if (ids.contains(rs.getInt(3))) {
                int arrayID = rs.getInt(1);
                
                if(cachedArrays.containsKey(arrayID)) { 
                    values.add(cachedArrays.get(arrayID));
                    
                } else {
                    ArrayDesign ad = new ArrayDesign(rs);
                    cacheArrayDesign(ad);
                    values.add(ad);
                }
            }
        }
        
        rs.close();
        s.close();
        return values;
    }
    
    public Collection<ArrayDesign> loadAllArrayDesigns() throws SQLException {
        HashSet<ArrayDesign> values = new HashSet<ArrayDesign>();
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id, name, genome from arraydesign");        
        
        while (rs.next()) {
            int arrayID = rs.getInt(1);
            ArrayDesign ad = null;
            if(cachedArrays.containsKey(arrayID)) { 
                ad = cachedArrays.get(arrayID);
            } else { 
                ad = new ArrayDesign(rs);
                cacheArrayDesign(ad);
            }
            values.add(ad);
        }
        
        rs.close();
        s.close();
        return values;
    }
    public GALFile loadGALFile(int dbid) throws SQLException {
        if(cachedGALFiles.containsKey(dbid)) { return cachedGALFiles.get(dbid); }
        GALFile gf = null;
        Statement s = cxn.createStatement();

        ResultSet rs = s.executeQuery("select id, name from galfiles where id=" + dbid);
        if(rs.next()) {
            gf = new GALFile(rs);
            cachedGALFiles.put(gf.getDBID(), gf);
            nameCachedGALFiles.put(gf.getName(), gf);
        }
        
        rs.close();
        s.close();
        
        if(gf == null) { throw new IllegalArgumentException("Unknown GALFile dbid: " + dbid); }
        return gf;        
    }
    
    public GALFile loadGALFile(String name) throws SQLException {
        name = name.replaceAll("^.*/","");
        if (name.matches("\\d+_D.*")) {
            name = name.replaceAll("_D.*","") + ".tdt";
        }
        if (name.matches("\\d+_2_D.*")) {
            name = name.replaceAll("_2_D.*","") + ".tdt";
        }

        if(nameCachedGALFiles.containsKey(name)) { return nameCachedGALFiles.get(name); }
        GALFile galfile = null;
        Statement s = cxn.createStatement();
        
        ResultSet rs = s.executeQuery("select id, name from galfiles where name='" + name + "'");
        if(rs.next()) { 
            galfile = new GALFile(rs);
            rs.close();
        } else { 
            rs.close();
            String seqName = "galfile_id";
            String idsql = Sequence.getInsertSQL(cxn, seqName);
            s.executeUpdate("insert into galfiles (id, name) values (" + idsql + ",'" + name + "')");
            
            rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, seqName));
            if(rs.next()) {
                int galID = rs.getInt(1);
                galfile = loadGALFile(galID);
            }
            rs.close();
        }

        s.close();
        
        if(galfile == null) { throw new IllegalStateException("Couldn't insert or file galfile \"" + name + "\""); }
        
        nameCachedGALFiles.put(galfile.getName(), galfile);
        cachedGALFiles.put(galfile.getDBID(), galfile);
        
        return galfile;
    }
    
    public ArrayDesign loadArrayDesign(int dbid) throws SQLException {
        if(cachedArrays.containsKey(dbid)) { return cachedArrays.get(dbid); }
        
        ArrayDesign ad = null;
        Statement s = cxn.createStatement();

        ResultSet rs = s.executeQuery("select id, name, genome from arraydesign where id=" + dbid);
        if(rs.next()) { ad = new ArrayDesign(rs); }
        
        rs.close();
        s.close();
        
        if(ad == null) { throw new IllegalArgumentException("Unknown ArrayDesign dbid: " + dbid); }
        
        cacheArrayDesign(ad);
        return ad;
    }
    
    public ArrayDesign loadArrayDesign(String name, Genome genome) throws SQLException {
        String seqName = "arraydesign_id";
        ArrayDesign result = retrieveArrayDesignFromCache(name, genome);
        if(result != null) { return result; }

        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery("select id, name, genome from arraydesign where name = '" + name + "'");        
        if (rs.next()) {
            try {
                result = new ArrayDesign(rs);
            } catch (UnknownRoleException ex) {
                throw new DatabaseException(ex.toString(),ex);
            }
        } else {
            rs.close();
            s.close();
            PreparedStatement ps = cxn.prepareStatement("insert into ArrayDesign (id,name,genome) values (" +
                                                        Sequence.getInsertSQL(cxn,seqName) + ",?,?)");
            ps.setString(1,name);
            ps.setInt(2,genome.getDBID());
            ps.execute();
            ps.close();
            s = cxn.createStatement();
            rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, seqName));
            if (rs.next()) {
                int arrayID = rs.getInt(1);
                result = loadArrayDesign(arrayID);
            }
        }
        rs.close();
        s.close();
        
        if(result == null) { 
            throw new IllegalStateException("Couldn't insert or find ArrayDesign " + name + "," + genome.getVersion()); 
        }
        
        cacheArrayDesign(result);
        return result;
    }
    public ArrayDesign getArrayDesign(String name, Genome genome) throws SQLException, NotFoundException {
        String seqName = "arraydesign_id";
        ArrayDesign result = retrieveArrayDesignFromCache(name, genome);

        if(result == null ) { 
        	Statement s = cxn.createStatement();
        	ResultSet rs = s.executeQuery("select id, name, genome from arraydesign where name = '" + name + "'");        
        	if (rs.next()) {
        		try {
        			result = new ArrayDesign(rs);
        		} catch (UnknownRoleException ex) {
        			throw new DatabaseException(ex.toString(),ex);
        		}
        	} else {
        		throw new NotFoundException("Couldn't find design " + name);
        	}
        }

        return result;
    }


    private FragDist parseFragDist(ResultSet rs) throws SQLException {
        FragDist dist = new FragDist();
        dist.dbid = rs.getInt(1);
        dist.name = rs.getString(2);
        dist.version = rs.getString(3);
        dist.description = rs.getString(4);
        PreparedStatement length = cxn.prepareStatement("select max(distance) from fragdistentry where distribution = ?");
        length.setInt(1,dist.dbid);
        ResultSet localrs = length.executeQuery();
        if (localrs.next()) {
            dist.values = new double[localrs.getInt(1) + 1];
            for (int i = 0; i < dist.values.length; i++) {
                dist.values[i] = 0;
            }
            localrs.close();
            PreparedStatement entries = cxn.prepareStatement("select distance,value from fragdistentry where distribution = ?");
            entries.setInt(1,dist.dbid);
            localrs = entries.executeQuery();
            while (localrs.next()) {
                dist.values[localrs.getInt(1)] = localrs.getDouble(2);
            }
        }else {            
            dist.values = new double[0];
        }
        localrs.close();
        length.close();
        return dist;

    }
    
    public Collection<FragDist> loadExperimentFragDists(String name, String version) 
    	throws SQLException, NotFoundException { 

    	LinkedList<FragDist> dists = new LinkedList<FragDist>();
    	String stmt = "select distinct(fragdist) from experiment where name=? and version=?"; 
    	PreparedStatement ps = cxn.prepareStatement(stmt);
    	
    	ps.setString(1, name);
    	ps.setString(2, version);
    	ResultSet rs = ps.executeQuery();
    	
    	while(rs.next()) { 
    		int id = rs.getInt(1);
    		FragDist fd = loadFragDist(id);
    		dists.addLast(fd);
    	}
    	
    	rs.close();
    	ps.close();
    	return dists;
    }
    
    /* retrieves the FragDist specified by the name and version.  Throws a NotFoundException if the Fragdist does not already
       exist in the database */
    public FragDist loadFragDist(String name, String version) throws SQLException, NotFoundException {
        PreparedStatement ps = cxn.prepareStatement("select id, name, version, description from fragdist where " + 
                                                    " name = ? and version = ?");
        ps.setString(1,name);
        ps.setString(2,version);
        ResultSet rs = ps.executeQuery();
        if (rs.next()) {
            FragDist result = parseFragDist(rs);
            rs.close();
            ps.close();
            return result;
        } else {
            rs.close();
            ps.close();
            throw new NotFoundException("Can't find fragdist " + name + ","+ version);
        }

    }
    public FragDist loadFragDist(int dbid) throws SQLException, NotFoundException {
        PreparedStatement ps = cxn.prepareStatement("select id, name, version, description from fragdist where " + 
                                                    " id = ?");
        ps.setInt(1,dbid);
        ResultSet rs = ps.executeQuery();
        if (rs.next()) {
            FragDist result = parseFragDist(rs);
            rs.close();
            ps.close();
            return result;
        } else {
            rs.close();
            ps.close();
            throw new NotFoundException("Can't find fragdist with id " + dbid);
        }
    }
    /* If the FragDist specified by the name and version already exists in the database, that FragDist is returned and
       description and values are ignored.  If it does not already exist, then name, version, description, and
       values are used to create a new database record */
    public FragDist getFragDist(String name, String version, String description, double values[]) throws SQLException {
        FragDist result;
        try {
            result = loadFragDist(name,version);
        } catch (NotFoundException e) {
            PreparedStatement ps = cxn.prepareStatement("insert into fragdist(id,name,version,description) values(" +
                                                        Sequence.getInsertSQL(cxn,"fragdist_id") + ",?,?,?)");
            ps.setString(1,name);
            ps.setString(2,version);
            ps.setString(3,description);
            ps.execute();
            ps.close();
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery(Sequence.getLastSQLStatement(cxn,"fragdist_id"));
            rs.next();
            result = new FragDist();
            result.dbid = rs.getInt(1);
            rs.close();
            stmt.close();
            result.name = name;
            result.version = version;
            result.description = description;
            ps = cxn.prepareStatement("insert into fragdistentry(distribution,distance,value) values(?,?,?)");
            ps.setInt(1,result.dbid);
            for (int i = 0; i < values.length; i++) {
                ps.setInt(2,i);
                ps.setDouble(3,values[i]);
                ps.execute();
            }
            result.values = values;
        }
        return result;

    }

    public Experiment loadExperiment(String name, String version, String replicate) throws SQLException, NotFoundException {
        PreparedStatement s = cxn.prepareStatement("select id, fragdist, species, factorone, cellsone, conditionone, " + 
                                                   "factortwo, cellstwo, conditiontwo, active from experiment where name = ? " + 
                                                   "and version = ? and replicate = ?");
        s.setString(1,name);
        s.setString(2,version);
        s.setString(3,replicate);
        ResultSet rs = s.executeQuery();
        Experiment result;
        if (rs.next()) {
            result = new Experiment();
            result.dbid = rs.getInt(1);
            result.name = name;
            result.version = version;
            result.replicate = replicate;
            result.fragdist = rs.getInt(2);
            result.species = rs.getInt(3);
            result.factorone = rs.getInt(4);
            result.cellsone = rs.getInt(5);
            result.conditionone = rs.getInt(6);
            result.factortwo = rs.getInt(7);
            result.cellstwo = rs.getInt(8);
            result.conditiontwo = rs.getInt(9);
            result.active = rs.getInt(10) > 0;
        } else {
            throw new NotFoundException("Couldn't find experiment with name = " + name + ", version = " + version + 
                                        " and replicate = " + replicate);
        }
        rs.close();
        s.close();
        return result;
    }

    public Experiment loadExperiment(int dbid) throws SQLException, NotFoundException {
        PreparedStatement s = cxn.prepareStatement("select name, version, replicate, fragdist, species, factorone, cellsone, conditionone, " + 
                                      "factortwo, cellstwo, conditiontwo, active from experiment where id = ?");
        s.setInt(1,dbid);
        ResultSet rs = s.executeQuery();
        Experiment result;
        if (rs.next()) {
            result = new Experiment();
            result.dbid = dbid;
            result.name = rs.getString(1);
            result.version = rs.getString(2);
            result.replicate = rs.getString(3);
            result.fragdist = rs.getInt(4);
            result.species = rs.getInt(5);
            result.factorone = rs.getInt(6);
            result.cellsone = rs.getInt(7);
            result.conditionone = rs.getInt(8);
            result.factortwo = rs.getInt(9);
            result.cellstwo = rs.getInt(10);
            result.conditiontwo = rs.getInt(11);
            result.active = rs.getInt(12) > 0;
        } else {
            throw new NotFoundException("Couldn't find experiment with id = " + dbid);
        }
        rs.close();
        s.close();
        return result;
    }

    public Collection<Experiment> loadExperiment(ExptNameVersion env) throws SQLException, NotFoundException {
        if (env.getReplicate() == null) {
            return loadExperiment(env.getName(),
                                  env.getVersion());
        } else {
            ArrayList<Experiment> output = new ArrayList<Experiment>();
            output.add(loadExperiment(env.getName(),
                                      env.getVersion(),
                                      env.getReplicate()));
            return output;
        }
    }

    public Collection<Experiment> loadExperiment(String name, String version)  throws SQLException, NotFoundException {
        PreparedStatement s = cxn.prepareStatement("select id, replicate, fragdist, species, factorone, cellsone, conditionone, " + 
                                                   "factortwo, cellstwo, conditiontwo, active from experiment where name = ? " + 
                                                   "and version = ?");
        s.setString(1,name);
        s.setString(2,version);
        ResultSet rs = s.executeQuery();
        ArrayList<Experiment> results = new ArrayList<Experiment>();
        while (rs.next()) {
            Experiment result = new Experiment();
            result.dbid = rs.getInt(1);
            result.name = name;
            result.version = version;
            result.replicate = rs.getString(2);
            result.fragdist = rs.getInt(3);
            result.species = rs.getInt(4);
            result.factorone = rs.getInt(5);
            result.cellsone = rs.getInt(6);
            result.conditionone = rs.getInt(7);
            result.factortwo = rs.getInt(8);
            result.cellstwo = rs.getInt(9);
            result.conditiontwo = rs.getInt(10);
            result.active = rs.getInt(11) > 0;
            results.add(result);
        }
        rs.close();
        s.close();
        return results;
    }

    /* throws NotFoundException if the species or fragdist cannot be found in the DB.
       DOES NOT throw a NotFoundException if the factors, cells, or conditions can't be found;
       they wil be created */
    public Experiment getExperiment(String name,
                                    String version,
                                    String replicate,
                                    String species,
                                    String fragdistname,
                                    String fragdistversion,
                                    String factorone,
                                    String factortwo, 
                                    String cellsone,
                                    String cellstwo,
                                    String conditionone,
                                    String conditiontwo,
                                    boolean active) throws SQLException, NotFoundException {
        Experiment expt;
        try {
            expt = loadExperiment(name,version,replicate);
        } catch (NotFoundException e) {
            expt = new Experiment();
            expt.name = name;
            expt.version = version;
            expt.replicate = replicate;
            Organism org = new Organism(species);
            expt.species = org.getDBID();
            MetadataLoader loader = new MetadataLoader();
            expt.fragdist = loadFragDist(fragdistname,fragdistversion).getDBID();
            expt.factorone = loader.getFactor(factorone).getDBID();
            expt.factortwo = loader.getFactor(factortwo).getDBID();
            expt.cellsone = loader.getCells(cellsone).getDBID();
            expt.cellstwo = loader.getCells(cellstwo).getDBID();
            expt.conditionone = loader.getCondition(conditionone).getDBID();
            expt.conditiontwo = loader.getCondition(conditiontwo).getDBID();
            expt.active = active;
            PreparedStatement ps = cxn.prepareStatement("insert into experiment(id,name,version,replicate,fragdist,species,cellsone,"+
                                                        "conditionone,factorone,cellstwo,conditiontwo,factortwo,active) values("+
                                                        Sequence.getInsertSQL(cxn,"experiment_id")+",?,?,?,?,?,?,?,?,?,?,?,?)");
            ps.setString(1,name);
            ps.setString(2,version);
            ps.setString(3,replicate);
            ps.setInt(4,expt.fragdist);
            ps.setInt(5,expt.species);
            ps.setInt(6,expt.cellsone);
            ps.setInt(7,expt.conditionone);
            ps.setInt(8,expt.factorone);
            ps.setInt(9,expt.cellstwo);
            ps.setInt(10,expt.conditiontwo);
            ps.setInt(11,expt.factortwo);
            ps.setInt(12,active ? 1 : 0);
            ps.execute();
            ps.close();
            Statement s = cxn.createStatement();
            ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn,"experiment_id"));
            rs.next();
            expt.dbid = rs.getInt(1);
            rs.close();
            s.close();
        }
        return expt;
    }

    public void mapExptToGenome(int exptid,
                                int genomeid) throws SQLException {
        PreparedStatement ps = cxn.prepareStatement("select count(*) from exptToGenome where experiment = ? and genome = ?");
        ps.setInt(1,exptid);
        ps.setInt(2,genomeid);
        ResultSet rs = ps.executeQuery();
        rs.next();
        if (rs.getInt(1) == 0) {
            rs.close();
            ps.close();
            ps = cxn.prepareStatement("insert into exptToGenome(experiment,genome) values(?,?)");
            ps.setInt(1,exptid);
            ps.setInt(2,genomeid);
            ps.execute();
        }
        rs.close();
        ps.close();
    }

    public JBDAnalysis getJBDAnalysis(int dbid) throws SQLException, NotFoundException {
        PreparedStatement ps = cxn.prepareStatement("select species, name, version, active from bayesanalysis where id = ?");
        ps.setInt(1,dbid);
        ResultSet rs = ps.executeQuery();
        JBDAnalysis analysis;
        if (rs.next()) {
            analysis = new JBDAnalysis(new Organism(rs.getInt(1)),
                                       rs.getString(2),
                                       rs.getString(3),
                                       dbid);
        } else {
            rs.close();
            ps.close();
            throw new NotFoundException("Couldn't find jbd analysis with id " + dbid);
        }
        rs.close();
        ps.close();
        ps = cxn.prepareStatement("select name, value from bayesparameters where analysis = ?");
        ps.setInt(1,dbid);
        rs = ps.executeQuery();
        HashMap<String,String> params = new HashMap<String,String>();
        while (rs.next()) {
            params.put(rs.getString(1),
                       rs.getString(2));
        }
        rs.close();
        ps.close();
        analysis.params = params;
        ps = cxn.prepareStatement("select experiment from bayesanalysisinputs where analysis = ?");
        ps.setInt(1,dbid);
        rs = ps.executeQuery();
        Set<Experiment> inputs = new HashSet<Experiment>();
        while (rs.next()) {
            inputs.add(loadExperiment(rs.getInt(1)));
        }
        analysis.inputs = inputs;        
        return analysis;
    }

    public void mapJBDAnalysisToGenome(int analysisid,
                                       int genomeid) throws SQLException {
        PreparedStatement ps = cxn.prepareStatement("select count(*) from bayesToGenome where analysis = ? and genome = ?");
        ps.setInt(1,analysisid);
        ps.setInt(2,genomeid);
        ResultSet rs = ps.executeQuery();
        rs.next();
        if (rs.getInt(1) == 0) {
            rs.close();
            ps.close();
            ps = cxn.prepareStatement("insert into bayesToGenome(analysis,genome) values(?,?)");
            ps.setInt(1,analysisid);
            ps.setInt(2,genomeid);
            ps.execute();
        }
        rs.close();
        ps.close();
    }

    public JBDAnalysis getJBDAnalysis(int species, String name, String version) throws SQLException, NotFoundException {
        PreparedStatement ps = cxn.prepareStatement("select id from bayesanalysis where species = ? and " +
                                                    " name = ?  and version = ?");
        ps.setInt(1,species);
        ps.setString(2,name);
        ps.setString(3,version);
        ResultSet rs = ps.executeQuery();
        if (rs.next()) {
            int dbid = rs.getInt(1);
            rs.close();
            ps.close();
            return getJBDAnalysis(dbid);
        } else {
            rs.close();
            ps.close();
            throw new NotFoundException("Couldn't find JBD analysis for " + species + "," + name + "," + version);
        }
    }

    public JBDAnalysis loadJBDAnalysis(int species, String name, String version) throws SQLException, NotFoundException {
        PreparedStatement ps = cxn.prepareStatement("select id from bayesanalysis where species = ? and " +
                                                    " name = ?  and version = ?");
        ps.setInt(1,species);
        ps.setString(2,name);
        ps.setString(3,version);
        ResultSet rs = ps.executeQuery();
        if (rs.next()) {
            int dbid = rs.getInt(1);
            rs.close();
            ps.close();
            return getJBDAnalysis(dbid);
        } else {
            rs.close();
            ps.close();
            ps = cxn.prepareStatement("insert into bayesanalysis(id,species,name,version,active) values  " + 
                                      "(" + Sequence.getInsertSQL(cxn,"analysis_id") + ",?,?,?,1)");
            ps.setInt(1,species);
            ps.setString(2,name);
            ps.setString(3,version);
            ps.execute();
            Statement s = cxn.createStatement();
            rs = s.executeQuery(Sequence.getLastSQLStatement(cxn,"analysis_id"));
            rs.next();
            JBDAnalysis analysis = new JBDAnalysis(new Organism(species),
                                                   name,
                                                   version,
                                                   rs.getInt(1));
            rs.close();
            s.close();
            ps.close();
            ps = cxn.prepareStatement("select name, value from bayesparameters where analysis = ?");
            ps.setInt(1,analysis.getDBID());
            rs = ps.executeQuery();
            Map<String,String> params = new HashMap<String,String>();
            while (rs.next()) {
                params.put(rs.getString(1), rs.getString(2));
            }
            analysis.params = params;
            rs.close();
            ps.close();

            ps = cxn.prepareStatement("select experiment from bayesanalysisinputs where analysis = ?");
            ps.setInt(1,analysis.getDBID());
            rs = ps.executeQuery();
            Set<Experiment> inputs = new HashSet<Experiment>();
            while (rs.next()) {
                inputs.add(loadExperiment(rs.getInt(1)));
            }
            analysis.inputs = inputs;


            return analysis;
        }
    }

    public void addJBDParam(JBDAnalysis analysis, String key, String value) throws SQLException {
        PreparedStatement ps = cxn.prepareStatement("insert into bayesparameters(analysis,name,value) values(?,?,?)");
        ps.setInt(1,analysis.getDBID());   
        ps.setString(2,key);
        ps.setString(3,value);
        ps.execute();
        ps.close();
    }
    public void addJBDParam(JBDAnalysis analysis, Map<String,String> params) throws SQLException {
        PreparedStatement ps = cxn.prepareStatement("insert into bayesparameters(analysis,name,value) values(?,?,?)");
        ps.setInt(1,analysis.getDBID());
        for (String k : params.keySet()) {
            if (analysis.params.containsKey(k)) {
                continue;
            }
            ps.setString(2,k);
            ps.setString(3,params.get(k));
            ps.execute();
            analysis.params.put(k,params.get(k));
        }
        ps.close();
    }
    public void addJBDInput(JBDAnalysis analysis, Experiment expt) throws SQLException {
        PreparedStatement ps = cxn.prepareStatement("insert into bayesanalysisinputs(analysis,experiment) values (?,?)");
        ps.setInt(1,analysis.getDBID());
        ps.setInt(2,expt.getDBID());
        ps.execute();
        ps.close();
    }
}
