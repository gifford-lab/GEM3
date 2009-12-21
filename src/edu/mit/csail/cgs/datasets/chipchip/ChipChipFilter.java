package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Cells;
import edu.mit.csail.cgs.datasets.general.Condition;
import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;

public class ChipChipFilter implements Closeable {
	
	private Genome genome;
	private java.sql.Connection cxn;
	
	public ChipChipFilter(Genome g) throws SQLException, UnknownRoleException {
		genome = g;
		cxn = DatabaseFactory.getConnection(ExptLocator.dbRole);
	}
	
	public Pair<Factor,Factor> getFactorData(MetadataLoader loader, ChipChipLocator loc) throws SQLException { 
		Pair<Factor,Factor> p = null;
		String query = "select e.factorone, e.factortwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			Factor cond1 = loader.loadFactor(c1);
			Factor cond2 = loader.loadFactor(c2);
			p = new Pair<Factor,Factor>(cond1, cond2);
		} 
		rs.close();
		ps.close();

		return p;		
	}
	
	public Pair<Condition,Condition> getConditionData(MetadataLoader loader, ChipChipLocator loc) throws SQLException { 
		Pair<Condition,Condition> p = null;
		String query = "select e.conditionone, e.conditiontwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			Condition cond1 = loader.loadCondition(c1);
			Condition cond2 = loader.loadCondition(c2);
			p = new Pair<Condition,Condition>(cond1, cond2);
		} else { 
		    System.err.println("Couldn't find any values for locator: " + loc.toString());
        }
		rs.close();
		ps.close();

		return p;		
	}
	
	public Pair<Cells,Cells> getCellsData(MetadataLoader loader, ChipChipLocator loc) throws SQLException {
		Pair<Cells,Cells> p = null;
		String query = "select e.cellsone, e.cellstwo from experiment e where e.name=? and e.version=?";
		PreparedStatement ps = cxn.prepareStatement(query);
		ps.setString(1, loc.getNameVersion().name);
		ps.setString(2, loc.getNameVersion().version);

		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			int c1 = rs.getInt(1), c2 = rs.getInt(2);
			Cells cells1 = loader.loadCells(c1);
			Cells cells2 = loader.loadCells(c2);
			p = new Pair<Cells,Cells>(cells1, cells2);
		} 
		rs.close();
		ps.close();

		return p;
	}

	public Collection<ExptLocator> 
		findBinding(Cells cells, Condition cond, Factor factor) throws SQLException {
		
		LinkedList<ExptLocator> locs = new LinkedList<ExptLocator>();
		Statement s = cxn.createStatement();
		
		String agilentQuery = "select e.name, e.version from experiment e, " +
				"exptToGenome e2g where e.id=e2g.experiment and e2g.genome=" + 
				genome.getDBID();
		
		if(cells != null) { 
			int dbid = cells.getDBID();
			agilentQuery += " and (e.cellsone=" + dbid + " or e.cellstwo=" + dbid + ")";
		}
		
		if(cond != null) { 
			int dbid = cond.getDBID();
			agilentQuery += " and (e.conditionone=" + dbid + " or e.conditiontwo=" + dbid + ")";
		}

		if(factor != null) { 
			int dbid = factor.getDBID();
			agilentQuery += " and (e.factorone=" + dbid + " or e.factortwo=" + dbid + ")";
		}
		
		ResultSet rs = s.executeQuery(agilentQuery);
		while(rs.next()) { 
			String n = rs.getString(1), v = rs.getString(2);
			ChipChipLocator loc = new ChipChipLocator(genome, n, v);
			locs.addLast(loc);
		}
		rs.close();
		
		s.close();
		return locs;
	}

	public void close() {
		DatabaseFactory.freeConnection(cxn);
		cxn=null;
	}

	public boolean isClosed() {
		return cxn==null;
	}
}
