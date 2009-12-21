package edu.mit.csail.cgs.datasets.alignments;

import java.sql.Clob;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import edu.mit.csail.cgs.utils.database.ClobHandler;
import edu.mit.csail.cgs.utils.database.Sequence;

public class InsertableAlignBlock { 
	
	public Integer id;
	public InsertableAlignment alignment;
	public Integer chromosome;
	public Integer start_pos, end_pos;
	public String strand;
	public char[] bit_string;
	public Integer gapped_length;
	
	public InsertableAlignBlock(InsertableAlignment align, int chrom, int start, int stop, char str, char[] bts, int gl) { 
		id = null;
		alignment = align;
		chromosome = chrom;
		start_pos = start;
		end_pos = stop;
		strand = String.valueOf(str);
		bit_string = bts;
		gapped_length = gl;
	}
	
	public static PreparedStatement prepareStatement(Connection cxn) throws SQLException { 
		String insert = String.format("insert into align_block " +
				"(id, alignment, chromosome, start_pos, end_pos, strand, bit_string, gapped_length) " +
				"values (%s, ?, ?, ?, ?, ?, ?, ?)",
				Sequence.getInsertSQL(cxn, "align_block_id"));
		return cxn.prepareStatement(insert);
	}

	public void setPreparedStatement(Connection cxn, PreparedStatement ps) throws SQLException { 
		ps.setInt(1, alignment.id);
		ps.setInt(2, chromosome);
		ps.setInt(3, start_pos);
		ps.setInt(4, end_pos);
		ps.setString(5, strand);
		ClobHandler.setClob(cxn, ps, 6, new String(bit_string));
		ps.setInt(7, gapped_length);
	}
	
	public int getInsertedID(Connection cxn, Statement s) throws SQLException { 
		ResultSet rs = s.executeQuery(Sequence.getLastSQLStatement(cxn, "align_block_id"));
		rs.next();
		id = rs.getInt(1);
		rs.close();
		return id;
	}
	
	public void updateClob(Statement s) throws SQLException { 
		ResultSet rs = s.executeQuery("select bit_string from align_block for update where id=" + id);
		Clob clob = rs.getClob(1);
		clob.setString((long)1, new String(bit_string));
		rs.close();
	}

	public boolean check() {
		System.out.println(String.format("chromosome=%d", chromosome));
		System.out.println(String.format("start_pos=%d", start_pos));
		System.out.println(String.format("end_pos=%d", end_pos));
		System.out.println(String.format("strand=%s", strand));
		System.out.println(String.format("gapped_length=%d", gapped_length));
		System.out.println(String.format("bit_string.length=%d", bit_string.length));
		return chromosome != null && start_pos != null && end_pos != null && end_pos >= start_pos && strand != null && bit_string != null && 
			bit_string.length > 0 && gapped_length != null && gapped_length == bit_string.length;
	}
}

