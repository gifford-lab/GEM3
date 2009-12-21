/**
 * 
 */
package edu.mit.csail.cgs.datasets.binding;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author Timothy Danford
 *
 */
public class BindingScanFilter implements Closeable {

	private MetadataLoader chipLoader;;
	private BindingScanLoader loader;
	private Genome genome;
	private boolean closed;
	
	public BindingScanFilter(Genome g, BindingScanLoader l) 
		throws SQLException, UnknownRoleException {
		
		genome = g;
		loader = l;
		chipLoader = new MetadataLoader();
		closed = false;
	}
	
	public boolean isClosed() { return closed; }
	
	public void close() { 
		chipLoader.close();
		loader = null;
		closed = true;
	}
	
	public Genome getGenome() { return genome; }
	public void setGenome(Genome g) { genome = g; }
	
	public Collection<Cells> getCells(BindingScan scan) 
		throws SQLException, UnknownRoleException {
		
		Collection<Pair<Integer,Integer>> exptIDs = loader.loadTypedExptPairs(scan);
		HashSet<Integer> dbids = new HashSet<Integer>();
		java.sql.Connection c = loader.getConnection();
		Statement s = c.createStatement();
		ResultSet rs;
		
		for(Pair<Integer,Integer> epair : exptIDs) { 
			int type = epair.getFirst(), eid = epair.getLast();

			switch(type) { 
			case BindingScan.AGILENT_TYPE:
				rs = s.executeQuery("select e.cellsone, e.cellstwo from experiment e, exptToGenome e2g" +
						" where e.id=e2g.experiment and e2g.genome=" + genome.getDBID() + 
						" and e.id=" + eid);
				while(rs.next()) { 
					int cellsone = rs.getInt(1), cellstwo = rs.getInt(2);
					dbids.add(cellsone);
					dbids.add(cellstwo);
				}
				break;
			case BindingScan.BAYES_TYPE:
				rs = s.executeQuery("select e.cellsone, e.cellstwo from experiment e, " +
						"bayesanalysisinputs a2e, bayesToGenome a2g where a2e.analysis=a2g.analysis" +
						" and a2e.analysis=" + eid + " and a2g.genome=" + genome.getDBID() + 
						" and a2e.experiment=e.id");
				while(rs.next()) { 
					int cellsone = rs.getInt(1), cellstwo = rs.getInt(2);
					dbids.add(cellsone);
					dbids.add(cellstwo);
				}
				break;
			case BindingScan.MSP_TYPE:
				rs = s.executeQuery("select e.cellsone, e.cellstwo from experiment e, " +
						"rosettaanalysisinputs a2e, rosettaToGenome a2g where a2e.analysis=a2g.analysis" +
						" and a2e.analysis=" + eid + " and a2g.genome=" + genome.getDBID() + 
						" and a2e.experiment=e.id");
				while(rs.next()) { 
					int cellsone = rs.getInt(1), cellstwo = rs.getInt(2);
					dbids.add(cellsone);
					dbids.add(cellstwo);
				}
				break;
			}
		}
		
		s.close();
		return chipLoader.loadAllCells(dbids);
	}
	
	public Collection<BindingScan> findScans(Cells cells, 
			Condition cond, Factor factor) throws SQLException  {
		
        if(genome==null) { return new LinkedList<BindingScan>(); }
        int gid = genome.getDBID();
        
        String tables1 = "bindingscanToExpt b2e, exptToGenome e2g";
        String tables2 = "bindingscanToExpt b2e, bayesToGenome a2g";
        String tables3 = "bindingscanToExpt b2e, rosettaToGenome a2g";
        
        String where1 = "b2e.scantype=" + BindingScan.AGILENT_TYPE + " and b2e.scanexpt=e2g.experiment and e2g.genome=" + gid;
        String where2 = "b2e.scantype=" + BindingScan.BAYES_TYPE + " and b2e.scanexpt=a2g.analysis and a2g.genome=" + gid;
        String where3 = "b2e.scantype=" + BindingScan.MSP_TYPE + " and b2e.scanexpt=a2g.analysis and a2g.genome=" + gid;
        
        if(cells != null || cond != null || factor != null) { 
            tables1 += ", experiment e";
            tables2 += ", bayesanalysisinputs a2e, experiment e";
            tables3 += ", rosettaanalysisinputs a2e, experiment e";
            
            where1 += " and b2e.scanexpt=e.id";
            where2 += " and b2e.scanexpt=a2e.analysis and a2e.experiment=e.id";
            where3 += " and b2e.scanexpt=a2e.analysis and a2e.experiment=e.id";
        }
        
        String query1 = "select b2e.scan from " + tables1 + " where " + where1;
        String query2 = "select b2e.scan from " + tables2 + " where " + where2;
        String query3 = "select b2e.scan from " + tables3 + " where " + where3;
        
		if(cells != null) { 
			int dbid = cells.getDBID();
			query1 += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
			query2 += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
			query3 += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
		}
		if(cond != null) { 
			int dbid = cond.getDBID();
			query1 += " and (e.conditionone=" + dbid + 
			" or e.conditiontwo=" + dbid + ")";
			query2 += " and (e.conditionone=" + dbid + 
			" or e.conditiontwo=" + dbid + ")";
			query3 += " and (e.conditionone=" + dbid + 
			" or e.conditiontwo=" + dbid + ")";
		}
		if(factor != null) { 
			int dbid = factor.getDBID();
			query1 += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
			query2 += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
			query3 += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
		}
		
		HashSet<Integer> dbids = new HashSet<Integer>();
		
		java.sql.Connection cxn = loader.getConnection();
		Statement s = cxn.createStatement();
		
        System.out.println("Querying agilent bindingscans...");
		ResultSet rs = s.executeQuery(query1);
		while(rs.next()) { 
			dbids.add(rs.getInt(1));
		}
		rs.close();
		
        System.out.println("Querying bayes bindingscans...");
		rs = s.executeQuery(query2);
		while(rs.next()) { 
			dbids.add(rs.getInt(1));
		}		
		rs.close();

        System.out.println("Querying rosetta bindingscans...");
		rs = s.executeQuery(query3);
		while(rs.next()) { 
			dbids.add(rs.getInt(1));
		}		
		rs.close();

		s.close();
        
        System.out.println("Done queryinq: loading scans by ID.");
		
		Collection<BindingScan> scans = loader.loadScans(genome, dbids);
		return scans;
	}
}
