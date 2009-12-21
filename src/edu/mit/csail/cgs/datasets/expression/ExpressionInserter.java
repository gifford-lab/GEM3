/*
 * Created on Mar 18, 2007
 */
package edu.mit.csail.cgs.datasets.expression;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.TimePoint;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class ExpressionInserter implements Closeable {
    
    public static String role = "expression";
    
    private java.sql.Connection cxn;

    private PreparedStatement insertExpt;
    private PreparedStatement insertExptParam;
    private PreparedStatement insertProbe;
    private PreparedStatement insertProbePlatform;
    private PreparedStatement insertMeasure;
    private PreparedStatement insertProcessing;
    private PreparedStatement insertProcessingParam;
    private PreparedStatement insertProcessingInput;
    private PreparedStatement insertSet, insertSetMember;

    public ExpressionInserter() throws SQLException {
        try {
            cxn = DatabaseFactory.getConnection(role);
        } catch (UnknownRoleException e) {
            throw new IllegalArgumentException("Unknown Role: " + role, e);
        }
        
        insertExpt = Experiment.prepareInsert(cxn);
        insertExptParam = Experiment.prepareInsertParam(cxn);
        insertProbe = Probe.prepareInsert(cxn);
        insertProbePlatform = ProbePlatform.prepareInsert(cxn);
        insertMeasure = ExprMeasurement.prepareInsert(cxn);
        insertProcessingParam = Processing.prepareInsertParam(cxn);
        insertProcessing = Processing.prepareInsert(cxn);
        insertProcessingInput = Processing.prepareInsertInput(cxn);
        insertSet = ExperimentSet.prepareInsert(cxn);
        insertSetMember = ExperimentSet.prepareInsertMember(cxn);
    }
    
    public ExpressionLoader getLoader() throws SQLException { return new ExpressionLoader(cxn); }

    public void close() {

        try {
            insertExpt.close();
            insertExptParam.close();
            insertProbe.close();
            insertProbePlatform.close();
            insertMeasure.close();
            insertProcessing.close();
            insertProcessingParam.close();
            insertProcessingInput.close();
            
            insertSet.close();
            insertSetMember.close();

        } catch (SQLException e) {
            e.printStackTrace();
        }

        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }
    
    public void beginTransaction() throws SQLException { 
        cxn.setAutoCommit(false);
    }
    
    public void commitTransaction() throws SQLException { 
        cxn.commit();
        cxn.setAutoCommit(true);
    }

    private int getLastID(String tableName) throws SQLException { 
        int id = -1;
        Statement s = cxn.createStatement();
        ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, tableName));
        if(rs.next()) { 
            id = rs.getInt(1);
        } else { 
            throw new IllegalStateException("Couldn't get last value of " + tableName + "?");
        }
        s.close();
        
        return id;              
    }
    
    private int getLastExperimentID() throws SQLException { return getLastID("experiment_id"); }
    private int getLastExperimentSetID() throws SQLException { return getLastID("experiment_set_id"); }
    private int getLastProbeID() throws SQLException { return getLastID("probe_id"); }
    private int getLastProbePlatformID() throws SQLException { return getLastID("probe_platform_id"); }
    private int getLastProcessingID() throws SQLException { return getLastID("processing_id"); }
    
    public int insertExperiment(String name, int type, boolean logscale, Cells cells, 
            Condition cond, TimePoint point, ProbePlatform platform) throws SQLException {
        insertExpt.setString(1, name);
        insertExpt.setInt(2, type);
        insertExpt.setInt(3, logscale ? 1 : 0);
        insertExpt.setInt(4, cells.getDBID());
        insertExpt.setInt(5, cond.getDBID());
        
        if(point != null) { 
            insertExpt.setInt(6, point.getDBID());
        } else { 
            insertExpt.setNull(6, java.sql.Types.INTEGER);
        }
        
        insertExpt.setInt(7, platform.getDBID());
        
        insertExpt.executeUpdate();
        return getLastExperimentID();
    }
    
    public void insertExperimentParam(Experiment expt, String k, String v) throws SQLException {
        if(expt.hasParam(k)) { 
            Statement s = cxn.createStatement();
            String update = "update expr_experiment_params set value='" + v + "' where key = '" + k + 
                "' and experiment=" + expt.getDBID();
            s.executeUpdate(update);
            s.close();
        } else { 
            insertExptParam.setInt(1, expt.getDBID());
            insertExptParam.setString(2, k);
            insertExptParam.setString(3, v);
            insertExptParam.executeUpdate();
        }
        expt.setParam(k, v);
    }
    
    public int insertExperimentSet(String name, String type) throws SQLException { 
        insertSet.setString(1, name);
        insertSet.setString(2, type);
        
        insertSet.executeUpdate();
        
        return getLastExperimentSetID();
    }
    
    public void insertExperimentSetMember(ExperimentSet set, Experiment member) throws SQLException {
        if(set.containsExperiment(member)) { throw new IllegalArgumentException("Already contains : " + member); }
        
        insertSetMember.setInt(1, set.getDBID());
        insertSetMember.setInt(2, member.getDBID());
        insertSetMember.executeUpdate();
        
        set.addMember(member);
    }

    public int insertProcessing(String type) throws SQLException { 
        insertProcessing.setString(1, type);
        
        insertProcessing.executeUpdate();
        
        return getLastProcessingID();
    }
    
    public void insertProcessingParam(Processing proc, String k, String v) throws SQLException {
        if(proc.hasParam(k)) { 
            Statement s = cxn.createStatement();
            String update = "update expr_processing_params set value='" + v + "' where key = '" + k + 
                "' and processing=" + proc.getDBID();
            s.executeUpdate(update);
            s.close();
        } else { 
            insertProcessingParam.setInt(1, proc.getDBID());
            insertProcessingParam.setString(2, k);
            insertProcessingParam.setString(3, v);
            insertProcessingParam.executeUpdate();
        }
        proc.setParam(k, v);
    }

    public void insertProcessingInputPair(Processing proc, Experiment input, Experiment output) 
    throws SQLException {
        if(proc.hasInputOutputPair(input, output)) { throw new IllegalArgumentException(); }

        insertProcessingInput.setInt(1, proc.getDBID());
        insertProcessingInput.setInt(2, input.getDBID());
        insertProcessingInput.setInt(3, output.getDBID());
        
        insertProcessingInput.executeUpdate();
        
        proc.addInputPair(input, output);
    }
    
    public int insertProbePlatform(String name, int type) throws SQLException { 
        insertProbePlatform.setString(1, name);
        insertProbePlatform.setInt(2, type);
        
        insertProbePlatform.executeUpdate();
        
        return getLastProbePlatformID();
    }
    
    public void insertProbePlatformGenomeMapping(ProbePlatform plat, Genome g) throws SQLException { 
    	PreparedStatement ps = cxn.prepareStatement("insert into " +
    			"probe_platform_to_genome " +
    			"(platform, genome) values (?, ?)");
    	
    	ps.setInt(1, plat.getDBID());
    	ps.setInt(2, g.getDBID());
    	ps.executeUpdate();
    	
    	ps.close();
    }
    
    public void insertProbes(Set<String> probeNames, int platformID) throws SQLException { 
        for(String name : probeNames) { 
            insertProbe.setString(1, name);
            insertProbe.setInt(2, platformID);
            insertProbe.executeUpdate();
        }
    }

	public void insertProbeLocations(Map<String,Set<Region>> probeLocations, 
			ProbePlatform platform) throws SQLException { 
		PreparedStatement pls = cxn.prepareStatement("insert into probe_location " + 
				"(probe, chromosome, startpos, stoppos) values (?, ?, ?, ?)");
		PreparedStatement fps = 
			cxn.prepareStatement("select id from probe where platform=? and name=?");

		for(String probeName : probeLocations.keySet()) {
			
			fps.setInt(1, platform.getDBID());
			fps.setString(2, probeName);
			ResultSet rs = fps.executeQuery();
			if(rs.next()) { 
				int probeID = rs.getInt(1);
				for(Region region : probeLocations.get(probeName)) { 
					pls.setInt(1, probeID);
					pls.setInt(2, region.getGenome().getChromID(region.getChrom()));
					pls.setInt(3, region.getStart());
					pls.setInt(4, region.getEnd());

					pls.executeUpdate();
				}
			} else { 
				System.err.println("Couldn't find probe \"" + probeName + "\"");
			}

			rs.close();
		}
		
		fps.close();
		pls.close();
	}
	
    public int insertProbe(String name, ProbePlatform platform) throws SQLException { 
        insertProbe.setString(1, name);
        insertProbe.setInt(2, platform.getDBID());
        
        insertProbe.executeUpdate();
        
        return getLastProbeID();
    }
    
    public void insertMeasurements(Experiment expt, Map<Probe,Double> values) throws SQLException { 
        for(Probe p : values.keySet()) { 
            insertMeasure.setInt(1, expt.getDBID());
            insertMeasure.setInt(2, p.getDBID());
            insertMeasure.setDouble(3, values.get(p));
            insertMeasure.executeUpdate();
        }
    }
    
    public void insertMeasurement(Experiment expt, Probe probe, Double value) throws SQLException { 
        insertMeasure.setInt(1, expt.getDBID());
        insertMeasure.setInt(2, probe.getDBID());
        insertMeasure.setDouble(3, value);
        
        insertMeasure.executeUpdate();
    }
}
