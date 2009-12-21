package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.io.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;

import edu.mit.csail.cgs.utils.database.*;

public class AlignmentInserter {
	
	public static void deleteAlignmentVersion(AlignmentVersion v) throws SQLException { 
		Connection cxn = DatabaseFactory.getConnection("alignment");
		boolean ac = cxn.getAutoCommit();
		cxn.setAutoCommit(false);
		
		Statement s = cxn.createStatement();
		
		int dbid = v.getDBID().getID();
		
		String deleteBlocks = String.format("delete from align_block where alignment in " +
					"(select id from alignment where version=%d)", dbid);
		String deleteAlignments = String.format("delete from alignment where version=%d", dbid);
		String deleteVersion = String.format("delete from alignment_version where id=%d", dbid);
		
		s.executeUpdate(deleteBlocks);
		s.executeUpdate(deleteAlignments);
		s.executeUpdate(deleteVersion);
		
		cxn.commit();

		s.close();
		cxn.setAutoCommit(ac);
	}

	public static void insertDataset(InsertableAlignmentDataset ds) throws SQLException { 
		Connection cxn = DatabaseFactory.getConnection("alignment");
		boolean ac = cxn.getAutoCommit();
		cxn.setAutoCommit(false);
		
		Statement s = cxn.createStatement();
		PreparedStatement vps = InsertableAlignmentVersion.prepareStatement(cxn);
		PreparedStatement aps = InsertableAlignment.prepareStatement(cxn);
		PreparedStatement bps = InsertableAlignBlock.prepareStatement(cxn);
		
		ds.version.setPreparedStatement(vps);
		vps.executeUpdate();
		ds.version.getInsertedID(cxn, s);
		
		System.out.println(String.format("Inserted version %s", ds.version.name));
		
		for(InsertableAlignment a : ds.alignments) { 
			a.setPreparedStatement(aps);
			aps.executeUpdate();
			a.getInsertedID(cxn, s);
			System.out.print("."); System.out.flush();
		}
		
		System.out.println(String.format("Inserted %d alignments ... ", ds.alignments.size()));
		
		for(InsertableAlignBlock b : ds.blocks) { 
			b.setPreparedStatement(cxn, bps);
			bps.executeUpdate();
			System.out.print("."); System.out.flush();
			b.getInsertedID(cxn, s);
		}
		
		System.out.println(String.format("Inserted %d blocks...", ds.blocks.size()));

		cxn.commit();
		cxn.setAutoCommit(ac);

		bps.close();
		aps.close();
		vps.close();
		s.close();
		
		DatabaseFactory.freeConnection(cxn);
	}
}
