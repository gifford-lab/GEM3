/*
 * Author: tdanford
 * Date: Nov 6, 2008
 */
package edu.mit.csail.cgs.tools.expression;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.models.*;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.database.*;

public class ProbePlatformTextLoader {
	
	public static void main(String[] args) {
		String name = Args.parseString(args, "name", null);
		String typename = Args.parseString(args, "type", null);
		String filename = Args.parseString(args, "probenames", null);
		
		if(name == null || typename == null || filename == null) {
			System.err.println("Usage: \n" +
					"ProbePlatformTextLoader --name=[platform-name] --type=[platform-type] --file=[input-file]\n" +
					"Where:\t[platform-name] is a string that hasn't been entered in the database already,\n" + 
					"      \t[platform-type] is PROBES if these probes are unique spots on an array,\n" +
					"      \t                or GENES if they are gene-summaries (and their names are gene names)\n" +
					"      \tand [input-file] is the path to a file which contains the list of probe names, one per line."); 
			return;
		}
		
		File file = new File(filename);
		int type = ProbePlatform.UNKNOWN_TYPE;
		if(typename.equals("GENES")) { type = ProbePlatform.GENE_TYPE; }
		if(typename.equals("PROBES")) { type = ProbePlatform.PROBE_TYPE; }
		
		if(type == ProbePlatform.UNKNOWN_TYPE) { 
			System.err.println(String.format("Error: Unknown probe type: \"%s\"", typename));
			return;
		}
		
		try {
			ProbePlatformTextLoader ttl = new ProbePlatformTextLoader(name, type, file);
			Connection cxn = DatabaseFactory.getConnection("expression");
			ttl.insert(cxn);
			DatabaseFactory.freeConnection(cxn);

		} catch (IOException e) {
			//e.printStackTrace();
			System.err.println(String.format("I/O Error: %s", e.getMessage()));
			
		} catch (UnknownRoleException e) {
			//e.printStackTrace();
			System.err.println(String.format("GSE Unknown Role: %s", e.getMessage()));
			System.err.println("(Have you correctly set up your expression_passwd file?)");
			
		} catch (SQLException e) {
			//e.printStackTrace();
			System.err.println(String.format("Database/Network Error: %s", e.getMessage()));
		}
	}
	
	private Map<String,ProbeModel> probes;
	private PlatformModel platform;
	
	public ProbePlatformTextLoader(String name, Integer type, File f) throws IOException { 
		DataFrame<ProbeModel> probemodels = 
			new DataFrame<ProbeModel>(ProbeModel.class, f, false, "name");
		probes = new HashMap<String,ProbeModel>();
		platform = new PlatformModel(name, type);
		
		for(int i = 0; i < probemodels.size(); i++) { 
			ProbeModel pm = probemodels.object(i);
			probes.put(pm.name, pm);
		}
		System.out.println(String.format(
				"Loaded %d probes from file %s", probes.size(), f.getName()));
	}
	
	public Collection<String> probeNames() { return probes.keySet(); }
	public String platformName() { return platform.name; }
	public Integer platformType() { return platform.type; }
	public Integer platformID() { return platform.id; }
	public Integer probeID(String pname) { return probes.get(pname).id; }
	
	public void insert(Connection cxn) throws SQLException {
		boolean ac = cxn.getAutoCommit();
		cxn.setAutoCommit(false);
		
		Statement s = cxn.createStatement();
		s.executeUpdate(String.format("insert into probe_platform (id, name, type) values " +
				"(%s, '%s', %d)", Sequence.getInsertSQL(cxn, "probe_platform_id"), 
				platform.name, platform.type));
		
		ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "probe_platform_id"));
		rs.next();
		platform.id = rs.getInt(1);
		rs.close();
		
		System.out.println(String.format("Inserted ProbePlatform -> %d", platform.id));
		
		PreparedStatement ps = cxn.prepareStatement("insert into probe (id, platform, name) " +
				"values (?, ?, ?)");

		int c = 0;
		for(String pname : probes.keySet()) { 
			ProbeModel model = probes.get(pname);
		
			ps.setString(1, Sequence.getInsertSQL(cxn, "probe_id"));
			ps.setInt(2, platform.id);
			ps.setString(3, model.name);
			
			ps.executeUpdate();
			
			rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "probe_id"));
			rs.next();
			model.id = rs.getInt(1);
			rs.close();
			
			c += 1;
			
			if(c % 1000 == 0) { 
				int k = c/1000;
				System.out.print(String.format(" %dk", k)); System.out.flush();
			}
		}
		
		System.out.println(String.format("\n# Probes Inserted: %d", c));
		
		ps.close();
		s.close();
		
		cxn.commit();
		cxn.setAutoCommit(ac);
		
		System.out.println("Finished.");
	}
	
	public static class PlatformModel extends Model {
		public Integer id;
		public String name;
		public Integer type;
		
		public PlatformModel(String n, Integer t) { 
			name = n; type = t;
		}
	}

	public static class ProbeModel extends Model {
		public Integer id;
		public String name;
	}
}
