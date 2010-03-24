package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.mit.csail.cgs.datasets.chippet.RunningOverlapSum;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class ChipSeqOverlapExpander implements Closeable, Mapper<StrandedRegion, RunningOverlapSum> {

    private ChipSeqLoader loader;
	private LinkedList<ChipSeqAlignment> alignments;
    private ChipSeqLocator locator;
    private java.sql.Connection cxn;
    private PreparedStatement stmt;
    private int extension;
    private Genome lastGenome;

    public ChipSeqOverlapExpander(ChipSeqLocator loc, int extension)  throws SQLException, IOException { 
        this.extension = extension;
        loader = new ChipSeqLoader();
        alignments = null;
        stmt = null;
        locator = loc;
    }   
    private void getAligns(Genome genome) throws SQLException {
        if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<ChipSeqAlignment>();
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
        cxn = DatabaseFactory.getConnection("chipseq");
        StringBuffer alignIDs = new StringBuffer();
        if (alignments.size() == 1) {
            alignIDs.append("alignment = " + alignments.get(0).getDBID());
        } else {

            alignIDs.append("alignment in (");
            for (int i = 0; i < alignments.size(); i++) {
                if (i == 0) {
                    alignIDs.append(alignments.get(i).getDBID());
                } else {
                    alignIDs.append("," + alignments.get(i).getDBID());
                }
            }
        }
        stmt = cxn.prepareStatement("select startpos, stoppos from chipseqhits where " + alignIDs.toString() +
                                    " and chromosome = ? and startpos > ? and stoppos < ? and strand = ?");
        stmt.setFetchSize(1000);
    }

	public RunningOverlapSum execute(StrandedRegion a) {
		try {
            Genome g = a.getGenome();
            getAligns(g);
            stmt.setInt(1, g.getChromID(a.getChrom()));
            stmt.setInt(2, a.getStart());
            stmt.setInt(3, a.getEnd());
            stmt.setString(4, a.getStrand() == '+' ? "+" : a.getStrand() == '-' ? "-" : " ");
            RunningOverlapSum sum = new RunningOverlapSum(g, a.getChrom());
            ResultSet rs = stmt.executeQuery();
            if (a.getStrand() == '+') {
                while (rs.next()) {
                    sum.addInterval(rs.getInt(1), rs.getInt(2) + extension);
                }
            } else if (a.getStrand() == '-') {
                while (rs.next()) {
                    sum.addInterval(rs.getInt(1) - extension, rs.getInt(2));
                }
            } 
            rs.close();
            return sum;
		} catch (SQLException e) {
			e.printStackTrace();
            throw new DatabaseException(e.toString(), e);
		}
	}
    public void close() {
		if(loader != null) { 
			loader.close();
            loader = null;
            if (alignments != null) { alignments.clear();}
            if (stmt != null) {
                try {
                    stmt.close();
                } catch (SQLException e) {
                    e.printStackTrace();
                } finally {
                    stmt = null;
                }
            }
            DatabaseFactory.freeConnection(cxn);
            cxn = null;
        }
	}
    public boolean isClosed() {
        return loader == null;
    }
}