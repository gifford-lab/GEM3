package edu.mit.csail.cgs.datasets.chipseq;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;

import java.sql.*;

/**
 * @author tdanford
 *
 * Represents a row from the chipseqalignments table, in the *chipseq schemas.
 * 
 * create sequence chipseqalignment_id;
 * create table chipseqalignments (
 *	id number(10),
 *	expt constraint fk_chipseqalignment_expt references chipseqexpts(id) not null,
 * 	name varchar2(1000) not null,
 * 	genome number(10) not null,
 *	constraint chipseqalignment_pk primary key (id)
 * );
 */
public class ChipSeqAlignment {
	
	private int dbid;
	private ChipSeqExpt expt;
	private String name;
	private Genome genome;
	
	public ChipSeqAlignment(ResultSet rs, ChipSeqLoader loader) throws SQLException, NotFoundException { 
		dbid = rs.getInt(1);
		int exptID = rs.getInt(2);
		name = rs.getString(3);
		int genomeID = rs.getInt(4);
		try {
			expt = loader.loadExperiment(exptID);
			genome = Organism.findGenome(genomeID);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
	
	ChipSeqAlignment(ResultSet rs, ChipSeqExpt ex) throws SQLException { 
		dbid = rs.getInt(1);
		expt = ex;  // exptID is still in position #2, we just don't need it.
		name = rs.getString(3);
		int genomeID = rs.getInt(4);
		try {
			genome = Organism.findGenome(genomeID);
		} catch (NotFoundException e) {
			e.printStackTrace();
		}		
	}
	/**
	 * This is a work around to use readdb id directly, without accessing the Oracle database
	 */
	ChipSeqAlignment(ChipSeqExpt expt, Genome genome){ 
		dbid = expt.getDBID();
		this.expt = expt;
		this.name = expt.getName();
		this.genome = genome;
	}

	public int getDBID() { return dbid; }
	public String getName() { return name; }
	public ChipSeqExpt getExpt() { return expt; }
	public Genome getGenome() { return genome; }
	
	public int hashCode() { 
		int code = 17;
		code += dbid; code *= 37;
		return code;
	}
	
	public String toString() { 
		return String.format("\"%s\" %s -> %s", name, expt.getName(), genome.getVersion());
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ChipSeqAlignment)) { return false; }
		ChipSeqAlignment a = (ChipSeqAlignment)o;
		if(dbid != a.dbid) { return false; }
		return true;
	}

	public static PreparedStatement createLoadByIDStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome from chipseqalignments where id=?";
		return c.prepareStatement(query);
	}
	
	public static PreparedStatement createLoadByNameAndExptStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome from chipseqalignments where name=? and expt=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByExptStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome from chipseqalignments where expt=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByGenomeStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome from chipseqalignments where genome=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createLoadAllByExptAndGenomeStatement(java.sql.Connection c) throws SQLException { 
		String query = "select id, expt, name, genome from chipseqalignments where expt=? and genome=?";
		return c.prepareStatement(query);
	}

	public static PreparedStatement createInsertStatement(java.sql.Connection c) throws SQLException { 
		String query = String.format(
				"insert into chipseqalignments (id, expt, name, genome) values (%s, ?, ?, ?)", 
				edu.mit.csail.cgs.utils.database.Sequence.getInsertSQL(c, "chipseqalignment_id"));
		return c.prepareStatement(query);
	}
	
}
