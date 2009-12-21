/*
 * Author: tdanford
 * Date: Oct 19, 2008
 */
package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ChipSeqReadGenerator implements Iterator<String>, Closeable {

	private Connection cxn;
	private Statement stmt;
	private ResultSet results;
	private String nextSeq;
	
	public ChipSeqReadGenerator(String expt, String repl) 
		throws SQLException, NotFoundException { 
		
		cxn = DatabaseFactory.getConnection("chipseq");
		PreparedStatement ps = 
			cxn.prepareStatement("select id from chipseqexpts where name=? and replicate=?");
		ps.setString(1, expt);
		ps.setString(2, repl);
		results = ps.executeQuery();
		int exptID = -1;
		if(results.next()) { 
			exptID = results.getInt(1);
		}
		results.close();
		ps.close();
		
		if(exptID == -1) { 
			throw new NotFoundException(String.format(
					"Unknown Expt/Replicate: \"%s\", \"%s\"", expt, repl));
		}
		
		stmt = cxn.createStatement();
		results = stmt.executeQuery(String.format(
				"select sequence from chipseqreads where expt=%d", exptID));
		
		findNextSequence();
	}
	
	private void findNextSequence() { 
		nextSeq = null;
		try {
			if(results.next()) { 
				nextSeq = results.getString(1);
			}
		} catch (SQLException e) {
			e.printStackTrace();
			nextSeq = null;
		}
		
		if(nextSeq == null && !isClosed()) { 
			close();
		}
	}

	public boolean hasNext() {
		if(isClosed()) { throw new IllegalStateException(); }
		
		return nextSeq != null;
	}

	public String next() {
		if(isClosed() || nextSeq == null) { 
			throw new NoSuchElementException();
		}
		
		String seq = nextSeq; 
		findNextSequence();
		return seq;
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}

	public void close() {
		try {
			results.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		DatabaseFactory.freeConnection(cxn);
		results = null;
		stmt = null;
		cxn = null;
		nextSeq = null;
	}

	public boolean isClosed() {
		return cxn == null;
	}
}
