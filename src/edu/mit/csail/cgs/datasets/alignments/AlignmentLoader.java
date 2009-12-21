/*
 * Author: tdanford
 * Date: Aug 22, 2008
 */
package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.io.*;
import java.sql.*;

import javax.sql.rowset.serial.SerialClob;

import edu.mit.csail.cgs.datasets.DBID;
import edu.mit.csail.cgs.datasets.SimpleDBID;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class AlignmentLoader implements edu.mit.csail.cgs.utils.Closeable {

	private java.sql.Connection alignCxn, coreCxn;
	private Statement alignStmt, coreStmt;
	
	private Map<Integer,String> cachedChroms;
	private Map<Integer,Genome> cachedChromGenomes;
	
	public AlignmentLoader() throws UnknownRoleException, SQLException { 
		alignCxn = DatabaseFactory.getConnection("alignment");
		coreCxn = DatabaseFactory.getConnection("core");
		
		alignStmt = alignCxn.createStatement();
		coreStmt = coreCxn.createStatement();
		
		cachedChroms = new HashMap<Integer,String>();
		cachedChromGenomes = new HashMap<Integer,Genome>();
	}
	
	public void close() {
		try { 
			alignStmt.close();
		} catch(SQLException e) {
			e.printStackTrace(System.err);
			alignStmt = null;
		}
		try {
			coreStmt.close();
		} catch (SQLException e) {
			e.printStackTrace(System.err);
			coreStmt = null;
		}
	}
	
	public boolean isClosed() { 
		return alignStmt == null || coreStmt == null;
	}
	
	private Genome lookupChromGenome(int chromID) throws SQLException { 
		if(!cachedChromGenomes.containsKey(chromID)) { 
			lookupChrom(chromID);
		}
		return cachedChromGenomes.get(chromID);
	}
	
	private String lookupChrom(int chromID) throws SQLException {
		if(cachedChroms.containsKey(chromID)) { 
			return cachedChroms.get(chromID);
		}
		String chrom = null;
		Genome g = null;
		String query = "select name, genome from chromosome where id=" + chromID;
		ResultSet rs = coreStmt.executeQuery(query);
		if(rs.next()) { 
			chrom = rs.getString(1);
			int genomeID = rs.getInt(2);
			
			try {
				g = Organism.findGenome(genomeID);
			} catch (NotFoundException e) {
				e.printStackTrace();
				throw new DatabaseException(e.getMessage());
			}
			
			cachedChroms.put(chromID, chrom);
			cachedChromGenomes.put(chromID, g);
		}
		rs.close();
		return chrom;
	}
	
	private Alignment loadAlignment(ResultSet rs) throws SQLException { 
		int alignID = rs.getInt(1);
		String params = rs.getString(2);
		int versionID = rs.getInt(3);
		Double score = rs.getDouble(4);
		
		return new CachedAlignment(alignID, versionID, params, score);
	}
	
	private Alignment loadAlignment(int id) throws SQLException { 
		Alignment value = null;		
		String query = "select id, params, version, score from alignment where id=" + id;

		ResultSet rs = alignStmt.executeQuery(query);
		if(rs.next()) {
			value = loadAlignment(rs);
		} 
		rs.close();
		
		return value;
	}
	
	private AlignBlock loadAlignBlock(ResultSet rs) throws SQLException {
		AlignBlock value = null;
		int blockID = rs.getInt(1);
		int alignID = rs.getInt(2);
		int chromID = rs.getInt(3);
		
		String chrom = lookupChrom(chromID);
		Genome genome = lookupChromGenome(chromID);
		
		int start = rs.getInt(4), end = rs.getInt(5);
		char strand = rs.getString(6).charAt(0);
		Clob bitclob = rs.getClob(7);
		int length = rs.getInt(8);
		
		Reader clobreader = bitclob.getCharacterStream();
		char[] bits = new char[length];
		int read = 0;
		try {
			read = clobreader.read(bits);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if(read >= length) { 
			value = new CachedAlignBlock(blockID, alignID, bits, 
					genome, chrom, start, end, strand, length);
		} else { 
			System.err.println(String.format("Read %d < gapped length %d", read, length));
		}		
		
		return value;
	}
	
	private AlignBlock loadAlignBlock(int id) throws SQLException {
		AlignBlock value = null;
		String query = "select id, alignment, chromosome, start_pos, " +
				"end_pos, strand, bit_string, gapped_length from align_block where id=" + id;
		ResultSet rs = alignStmt.executeQuery(query);
		if(rs.next()) {
			value = loadAlignBlock(rs);
		} 
		rs.close();
		
		return value;
	}

	private Collection<AlignBlock> loadAlignBlocks(Alignment a) throws SQLException {
		LinkedList<AlignBlock> blocks = new LinkedList<AlignBlock>();
		String query = "select ab.id, ab.alignment, ab.chromosome, ab.start_pos, " +
		"ab.end_pos, ab.strand, ab.bit_string, ab.gapped_length " +
		"from align_block ab where ab.alignment=?";
	
		PreparedStatement ps = alignCxn.prepareStatement(query);
		int aid = a.getDBID().getID();
		ps.setInt(1, aid);
		
		ResultSet rs = ps.executeQuery();
		while(rs.next()) { 
			AlignBlock block = loadAlignBlock(rs);
			if(block != null) { 
				blocks.add(block);
			}
		}
		rs.close();
		ps.close();
		
		return blocks;
	}
	
	private Collection<AlignBlock> loadAlignBlocks(AlignmentVersion version, Region r) throws SQLException { 
		Genome g = r.getGenome();
		int chromID = g.getChromID(r.getChrom());
		return loadAlignBlocks(version.getDBID(), chromID, r.getStart(), r.getEnd());
	}
	
	private Collection<AlignBlock> loadAlignBlocks(DBID versionID, int chromID, int start, int end) 
		throws SQLException {
		
		LinkedList<AlignBlock> blocks = new LinkedList<AlignBlock>();
		String query = "select ab.id, ab.alignment, ab.chromosome, ab.start_pos, " +
			"ab.end_pos, ab.strand, ab.bit_string, ab.gapped_length " +
			"from align_block ab, alignment where ab.alignment=alignment.id and " +
			"alignment.version=? and chromosome=? and " +
			"((start_pos >= ? and start_pos <= ?) or (start_pos <= ? and end_pos >= ?))";
		PreparedStatement ps = alignCxn.prepareStatement(query);
		
		int vid = versionID.getID();
		ps.setInt(1, vid);
		ps.setInt(2, chromID);
		ps.setInt(3, start);
		ps.setInt(4, end);
		ps.setInt(5, start);
		ps.setInt(6, start);
		
		ResultSet rs = ps.executeQuery();
		while(rs.next()) { 
			AlignBlock value = loadAlignBlock(rs);
			if(value != null) { blocks.addLast(value); }
		}
		rs.close();
		
		ps.close();
		return blocks;
	}
	
	private AlignmentVersion loadAlignmentVersion(ResultSet rs) throws SQLException { 
		int alignVersionID = rs.getInt(1);
		String name = rs.getString(2);
		return new CachedAlignmentVersion(alignVersionID, name);
	}
	
	private AlignmentVersion loadAlignmentVersion(int id) throws SQLException {
		AlignmentVersion value = null;
		String query = "select id, name from alignment_version where id=" + id;
		
		ResultSet rs = alignStmt.executeQuery(query);
		if(rs.next()) { 
			value = loadAlignmentVersion(rs);
		}
		rs.close();
		
		return value;
	}
	
	public Collection<String> getAlignmentVersionNames() throws SQLException { 
		LinkedList<String> names = new LinkedList<String>();
		ResultSet rs = alignStmt.executeQuery("select name from alignment_version order by name");
		while(rs.next()) { 
			names.addLast(rs.getString(1));
		}
		rs.close();
		return names;
	}
	
	public AlignmentVersion loadAlignmentVersion(String n) throws SQLException {
		AlignmentVersion value = null;
		String query = "select id, name from alignment_version where name=?";

		PreparedStatement ps = alignCxn.prepareStatement(query);
		ps.setString(1, n);
		
		ResultSet rs = ps.executeQuery();
		if(rs.next()) { 
			value = loadAlignmentVersion(rs);
		}

		rs.close();
		ps.close();
		return value;
	}
	
	private class CachedAlignBlock implements AlignBlock {
		
		private SimpleDBID dbid; 
		private int alignID;
		private char[] bitstring;
		private Genome genome;
		private String chrom;
		private int startpos, stoppos;
		private char strand;
		private int gappedLength;
		
		public CachedAlignBlock(
				int id, int alignid, char[] bits, Genome g, String chr, 
				int start, int stop, char str, int gap) {

			dbid = new SimpleDBID(id);
			alignID = alignid;
			bitstring = bits;
			genome = g;
			chrom = chr;
			startpos = start;
			stoppos = stop;
			strand = str;
			gappedLength = gap;
		}

		public Alignment getAlignment() {
			try {
				return loadAlignment(alignID);
			} catch (SQLException e) {
				e.printStackTrace();
				throw new DatabaseException(e.getMessage());
			}
		}

		public char[] getBitString() { return bitstring; }
		public String getChrom() { return chrom; }
		public DBID getDBID() { return dbid; }
		public int getGappedLength() { return gappedLength; }
		public int getStartPos() { return startpos; }
		public int getStopPos() { return stoppos; }
		public char getStrand() { return strand; }
		public Genome getGenome() { return genome; }
		
		public int hashCode() { return dbid.hashCode(); }
		
		public boolean equals(Object o) { 
			if(!(o instanceof CachedAlignBlock)) { return false; }
			return ((CachedAlignBlock)o).dbid.equals(dbid);
		}
	}
	
	private class CachedAlignment implements Alignment {
		
		private SimpleDBID dbid; 
		private int versionID;
		private String params; 
		private Double score;
		
		public CachedAlignment(int id, int vid, String p, Double sc) { 
			dbid = new SimpleDBID(id);
			versionID = vid;
			params = p;
			score = sc;
		}

		public DBID getDBID() { return dbid; }
		public String getParams() { return params; }
		public Double getScore() { return score; }

		public AlignmentVersion getVersion() {
			try {
				return loadAlignmentVersion(versionID);
			} catch (SQLException e) {
				e.printStackTrace();
				throw new DatabaseException(e.getMessage());
			}
		}

		public Collection<AlignBlock> getAlignBlocks() {
			try {
				return loadAlignBlocks(this);
			} catch (SQLException e) {
				e.printStackTrace();
				throw new DatabaseException(e.getMessage());
			}
		} 

		public int hashCode() { return dbid.hashCode(); }
		
		public boolean equals(Object o) { 
			if(!(o instanceof CachedAlignment)) { return false; }
			return ((CachedAlignment)o).dbid.equals(dbid);
		}
	}
	
	private class CachedAlignmentVersion implements AlignmentVersion {
		
		private SimpleDBID dbid; 
		private String name;
		
		public CachedAlignmentVersion(int id, String n) { 
			dbid = new SimpleDBID(id);
			name = n;
		}

		public DBID getDBID() { return dbid; }
		public String getName() { return name; }

		public Collection<AlignBlock> getAlignBlocks(Region r) {
			try {
				return loadAlignBlocks(this, r);
			} catch (SQLException e) {
				e.printStackTrace();
				throw new DatabaseException(e.getMessage());
			}
		}

		public int hashCode() { return dbid.hashCode(); }
		
		public boolean equals(Object o) { 
			if(!(o instanceof CachedAlignmentVersion)) { return false; }
			return ((CachedAlignmentVersion)o).dbid.equals(dbid);
		}
	}
}
