package edu.mit.csail.cgs.datasets.chipchip;

import java.io.*;
import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;

public class ChipChipDataset {

    public static int CHIP = 1, EXPRESSION = 2, CGH = 3, RULER = 4;
    
    public static List<String> getAvailableVersions(int genomeid) { 
        LinkedList<String> lst = new LinkedList<String>();
        java.sql.Connection cxn;
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement s = cxn.createStatement();
            String queryString = "select unique(e.version) from experiment e, arraydesign ad" +
                    ", genome g where e.design=ad.id and ad.genome=" + genomeid;
            ResultSet rs = s.executeQuery(queryString);
            while(rs.next()) { lst.addLast(rs.getString(1)); }
            rs.close();
            s.close();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            throw new IllegalArgumentException(se);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");
        }

        DatabaseFactory.freeConnection(cxn);
        return lst;
    }
    
    private Genome genome;
    private int speciesid;

    public ChipChipDataset (Genome g) {
        genome = g;
        java.sql.Connection cxn = null;
        Statement stmt = null;
        ResultSet rs = null;
        try {
            cxn = DatabaseFactory.getConnection("core");
            stmt = cxn.createStatement();
            rs = stmt.executeQuery("select species from genome where id = " + genome.getDBID());
            if (rs.next()) {
                speciesid = rs.getInt(1);
            } else {
                throw new DatabaseException("Couldn't find " + genome.getName());
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + genome.getName() + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core");
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }         
    }

    public List<ExptNameVersion> getExpts() {
        java.sql.Connection cxn = null;
        Statement stmt = null;
        ResultSet rs = null;
        List<ExptNameVersion> expts = new ArrayList<ExptNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            stmt = cxn.createStatement();
            rs = stmt.executeQuery("select distinct(name,version) from experiment where species = " + speciesid);
            while (rs.next()) {
                expts.add(new ExptNameVersion(rs.getString(1),rs.getString(2)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return expts;
    }
    public List<ExptNameVersion> getExptsReplicates() {
        java.sql.Connection cxn = null;
        Statement stmt = null;
        ResultSet rs = null;
        List<ExptNameVersion> expts = new ArrayList<ExptNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            stmt = cxn.createStatement();
            rs = stmt.executeQuery("select name,version,replicate from experiment where active = 1 and species = " + speciesid);
            while (rs.next()) {
                expts.add(new ExptNameVersion(rs.getString(1),rs.getString(2),rs.getString(3)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return expts;
    }
    public ChipChipData getData(String exptname, String version) throws NotFoundException {
        return new edu.mit.csail.cgs.datasets.chipchip.SQLData(exptname,version,genome.getDBID(),null);
    }
    public ChipChipData getData(ExptNameVersion env) throws NotFoundException {
        Set<String> repls;
        if (env.replicate != null) {
            repls = new HashSet<String>();
            repls.add(env.replicate);
        } else {
            repls = null;
        }
        return new edu.mit.csail.cgs.datasets.chipchip.SQLData(env.name,env.version,genome.getDBID(),repls);
    }
    public ChipChipData getWCE(String exptname, String version) throws NotFoundException {
        return new edu.mit.csail.cgs.datasets.chipchip.SQLDataWCEVals(exptname,version,genome.getDBID(),null);
    }
    public ChipChipData getWCE(ExptNameVersion env) throws NotFoundException {
        Set<String> repls;
        if (env.replicate != null) {
            repls = new HashSet<String>();
            repls.add(env.replicate);
        } else {
            repls = null;
        }
        return new edu.mit.csail.cgs.datasets.chipchip.SQLDataWCEVals(env.name,env.version,genome.getDBID(),repls);
    }
    public ChipChipCoeffs getCoeffs(String exptname, String version) throws NotFoundException {
        return new SQLCoeffs(exptname,version,speciesid);
    }

    public List<AnalysisNameVersion> getMLEAnalyses() {
        java.sql.Connection cxn = null;
        List<AnalysisNameVersion> mleanalyses = new ArrayList<AnalysisNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            String query = "select name, version from mleanalysis where active = 1 and " +
                "species = " + speciesid;
            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
                mleanalyses.add(new AnalysisNameVersion(rs.getString(1),rs.getString(2)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return mleanalyses;
    }
    /* get all of the MLE analyses derived from a particular experiment */
    public List<AnalysisNameVersion> getMLEAnalyses(NameVersion nv) {
        java.sql.Connection cxn = null;
        List<AnalysisNameVersion> mleanalyses = new ArrayList<AnalysisNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            String query = "select mleanalysis.name, mleanalysis.version from mleanalysis, mleanalysisinputs, experiment where " +
                "experiment.species = " + speciesid + " and experiment.version = '" +
                nv.version + "' and experiment.name = '" + nv.name + "' and experiment.id = mleanalysisinputs.experiment " +
                " and mleanalysisinputs.analysis = mleanalysis.id mleanalysis.active = 1";
            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
                mleanalyses.add(new AnalysisNameVersion(rs.getString(1),rs.getString(2)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return mleanalyses;
    }
    public ChipChipMLE getMLE(String analysis, String mleversion) throws NotFoundException {
        return new SQLMLE(analysis,mleversion,genome.getDBID());
    }
    public ChipChipMLE getMLE(AnalysisNameVersion anv) throws NotFoundException {
        return getMLE(anv.name,anv.version);
    }

    public Map<String,String> getMLEParams(String analysis, String resultversion) throws NotFoundException {
        java.sql.Connection cxn = null;
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select id from mleanalysis where " +
                                             " species = " + speciesid + 
                                             " and version = '" + resultversion + "'" + 
                                             " and name = '" + analysis + "'");
            rs.next();
            int id = rs.getInt(1);
            rs.close();
            rs = stmt.executeQuery("select name, value from mleparameters where analysis = " + id);
            Map<String,String> map = new HashMap<String,String>();
            while (rs.next()) {
                map.put(rs.getString(1),rs.getString(2));
            }
            rs.close();
            stmt.close();
            return map;
            
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }   
    }
    public List<ExptNameVersion> getMLEInputs(AnalysisNameVersion anv) {
        java.sql.Connection cxn = null;
        ArrayList<ExptNameVersion> envs = new ArrayList<ExptNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select id from mleanalysis where " +
                                             " species = " + speciesid + 
                                             " and version = '" + anv.version + "'" + 
                                             " and name = '" + anv.name + "'");
            rs.next();
            int id = rs.getInt(1);
            rs.close();
            rs = stmt.executeQuery("select experiment.name, experiment.version, experiment.replicate from " +
                                   "experiment, mleanalysisinputs where mleanalysisinputs.analysis = " + id +
                                   " and experiment.id = mleanlaysisinput.experiment");
            while (rs.next()) {
                envs.add(new ExptNameVersion(rs.getString(1),
                                             rs.getString(2),
                                             rs.getString(3)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }   
        return envs;
    }

    public List<AnalysisNameVersion> getBayesAnalyses() {
        java.sql.Connection cxn = null;
        List<AnalysisNameVersion> bayesanalyses = new ArrayList<AnalysisNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            String query = "select name, version from bayesanalysis where active = 1 and " +
                "experiment.species = " + speciesid;
            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
                bayesanalyses.add(new AnalysisNameVersion(rs.getString(1),rs.getString(2)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return bayesanalyses;
    }
    public List<AnalysisNameVersion> getBayesAnalyses(NameVersion nv) {
        java.sql.Connection cxn = null;
        List<AnalysisNameVersion> bayesanalyses = new ArrayList<AnalysisNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            String query = "select bayesanalysis.name, bayesanalysis.version from bayesanalysis, bayesanalysisinputs, experiment where " +
                "experiment.species = " + speciesid + " and experiment.version = '" +
                nv.version + "' and experiment.name = '" + nv.name + "' and experiment.id = bayesanalysisinputs.experiment " +
                " and bayesanalysisinputs.analysis = bayesanalysis.id and bayesanalysis.active = 1";
            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
                bayesanalyses.add(new AnalysisNameVersion(rs.getString(1),rs.getString(2)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        } 
        return bayesanalyses;
    }


    public ChipChipBayes getBayes(String analysis, String bayesversion) throws NotFoundException {
        return new SQLBayes(analysis,bayesversion,genome.getDBID());
    }
    public ChipChipBayes getBayes(AnalysisNameVersion anv) throws NotFoundException {
        return getBayes(anv.name,anv.version);
    }
    public Map<String,String> getBayesParams(String analysis, String resultversion) throws NotFoundException {
        java.sql.Connection cxn = null;
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select id from bayesanalysis where " +
                                             " species = " + speciesid + 
                                             " and version = '" + resultversion + "'" + 
                                             " and and name = '" + analysis + "'");
            rs.next();                
            int id = rs.getInt(1);
            rs.close();
            rs = stmt.executeQuery("select name, value from bayesparameters where analysis = " + id);
            Map<String,String> map = new HashMap<String,String>();
            while (rs.next()) {
                map.put(rs.getString(1),rs.getString(2));
            }
            rs.close();
            stmt.close();
            return map;
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }   
    }
    public List<ExptNameVersion> getBayesInputs(AnalysisNameVersion anv) {
        java.sql.Connection cxn = null;
        ArrayList<ExptNameVersion> envs = new ArrayList<ExptNameVersion>();
        try {
            cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select id from bayesanalysis where " +
                                             " species = " + speciesid + 
                                             " and version = '" + anv.version + "'" + 
                                             " and name = '" + anv.name + "'");
            rs.next();
            int id = rs.getInt(1);
            rs.close();
            rs = stmt.executeQuery("select experiment.name, experiment.version, experiment.replicate from " +
                                   "experiment, bayesalysisinputs where bayesanalysisinputs.analysis = " + id +
                                   " and experiment.id = bayesanlaysisinput.experiment");
            while (rs.next()) {
                envs.add(new ExptNameVersion(rs.getString(1),
                                             rs.getString(2),
                                             rs.getString(3)));
            }
            rs.close();
            stmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't find " + speciesid + ": "+ ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        } finally {
            DatabaseFactory.freeConnection(cxn);
        }   
        return envs;
    }

    public List<AnalysisNameVersion> getMSPAnalyses() { 
        List<AnalysisNameVersion>mspanalyses = new ArrayList<AnalysisNameVersion>();
        try { 
            java.sql.Connection c = DatabaseFactory.getConnection("chipchip");
            PreparedStatement s = c.prepareStatement("select name, version from rosettaanalysis where species=? and active = 1");
            s.setInt(1,speciesid);
            ResultSet rs = s.executeQuery();
            while(rs.next()) { 
                String expt = rs.getString(1);
                String version = rs.getString(2);
                AnalysisNameVersion env = new AnalysisNameVersion(expt, version);
                mspanalyses.add(env);
            }
            rs.close();
            s.close();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip");      
        }
        return mspanalyses;
    }
    public ChipChipMSP getMSP(String analysis, String version) throws NotFoundException {
        return new SQLMSP(analysis,version,genome.getDBID());
    }
    public ChipChipMSP getMSP(AnalysisNameVersion anv) throws NotFoundException {
        return getMSP(anv.name,anv.version);
    }

    public String getName () {
        return "ChipChip Data for " + genome.getName();
    }
    public Genome getGenome() { 
    	return genome;
    }
    public int getExptType(ExptNameVersion env) throws SQLException, NotFoundException {
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        MetadataLoader core = new MetadataLoader();
        Collection<Experiment> expts = loader.loadExperiment(env);
        int type = CHIP;
        for (Experiment expt : expts) {
            try {
                if (core.loadFactor(expt.getFactorOne()).getName().toLowerCase().matches(".*ruler.*")) {
                    type = RULER;
                }
                if (core.loadFactor(expt.getFactorOne()).getName().toLowerCase().matches(".*rna.*")) {
                    type = EXPRESSION;
                }
                if (core.loadFactor(expt.getFactorOne()).getName().toLowerCase().matches("cgh")) {
                    type = CGH;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        return type;
    }


}
