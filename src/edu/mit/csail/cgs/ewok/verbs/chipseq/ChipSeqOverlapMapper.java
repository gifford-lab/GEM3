package edu.mit.csail.cgs.ewok.verbs.chipseq;

import java.sql.*;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.chippet.RunningOverlapSum;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

public class ChipSeqOverlapMapper implements Closeable,  Mapper<Region, RunningOverlapSum> {

    private ChipSeqLoader loader;
	private LinkedList<ChipSeqAlignment> alignments;
    private java.sql.Connection cxn;
    private PreparedStatement stmtStranded, stmtBoth;
    private int extension;
    private int shift=0;
    private boolean shifting=false;

    public ChipSeqOverlapMapper(ChipSeqLocator loc, int extension) throws SQLException, IOException { 
    	this.extension = extension;
    	loader = new ChipSeqLoader();
    	alignments = new LinkedList<ChipSeqAlignment>();

    	try { 
    		alignments.addAll(loc.loadAlignments(loader));
    	} catch(SQLException e) { 
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
            alignIDs.append(")");
        }
        stmtStranded = cxn.prepareStatement("select startpos, stoppos from chipseqhits where " + alignIDs.toString() +
                                            " and chromosome = ? and startpos > ? and stoppos < ? and strand = ?");
        stmtBoth = cxn.prepareStatement("select startpos, stoppos, strand from chipseqhits where " + alignIDs.toString() +
                                        " and chromosome = ? and startpos > ? and stoppos < ?");
        stmtStranded.setFetchSize(5000);
        stmtBoth.setFetchSize(5000);
    }   

	public RunningOverlapSum execute(StrandedRegion a) {
		try {
            Genome g = a.getGenome();
            stmtStranded.setInt(1, g.getChromID(a.getChrom()));
            stmtStranded.setInt(2, a.getStart());
            stmtStranded.setInt(3, a.getEnd());
            stmtStranded.setString(4, a.getStrand() == '+' ? "+" : a.getStrand() == '-' ? "-" : " ");
            RunningOverlapSum sum = new RunningOverlapSum(g, a.getChrom());
            ResultSet rs = stmtStranded.executeQuery();
            if (a.getStrand() == '+') {
                while (rs.next()) {
                	if(shifting)
                		sum.addInterval(rs.getInt(1)+shift-(extension/2), rs.getInt(2) +shift+(extension/2));
                	else
                		sum.addInterval(rs.getInt(1), rs.getInt(2)+ extension);
                }
            } else if (a.getStrand() == '-') {
                while (rs.next()) {
                	if(shifting)
                		sum.addInterval(rs.getInt(1)-shift - (extension/2), rs.getInt(2)-shift+(extension/2));
                	else
                		sum.addInterval(rs.getInt(1)- extension, rs.getInt(2));
                }
            } 
            rs.close();
            return sum;
		} catch (SQLException e) {
			e.printStackTrace();
            throw new DatabaseException(e.toString(), e);
		}
	}
	public RunningOverlapSum execute(Region a) {
        if (a instanceof StrandedRegion) {
            return execute((StrandedRegion)a);
        }

		try {
            Genome g = a.getGenome();            
            stmtBoth.setInt(1, g.getChromID(a.getChrom()));
            stmtBoth.setInt(2, a.getStart());
            stmtBoth.setInt(3, a.getEnd());
            RunningOverlapSum sum = new RunningOverlapSum(g, a.getChrom());
            ResultSet rs = stmtBoth.executeQuery();
            while (rs.next()) {
                if (rs.getString(3).charAt(0) == '+') {
                	if(shifting)
                		sum.addInterval(rs.getInt(1)+shift-(extension/2), rs.getInt(2)+shift+(extension/2));
                	else
                		sum.addInterval(rs.getInt(1), rs.getInt(2)+ extension);
                } else {
                	if(shifting)
                		sum.addInterval(rs.getInt(1)-shift-(extension/2), rs.getInt(2)-shift+(extension/2));
                	else
                		sum.addInterval(rs.getInt(1)-extension, rs.getInt(2));
                }
            }
            rs.close();
            return sum;
		} catch (SQLException e) {
			e.printStackTrace();
            throw new DatabaseException(e.toString(), e);
		}
	}    
    public void setExtension (int e) {
        extension = e;
    }
    public void setShift (int s) {
        shift = s;
        if(s>0)
        	shifting=true;
        else
        	shifting=false;
    }
    public void close() {
		if(loader != null) { 
			loader.close();
            loader = null;
            alignments.clear();
            try {
                stmtStranded.close();
                stmtBoth.close();
            } catch (SQLException e) {
                e.printStackTrace();
            } finally {
                stmtStranded = null;
                stmtBoth = null;
            }
            DatabaseFactory.freeConnection(cxn);
            cxn = null;
        }
	}
    public boolean isClosed() {
        return loader == null;
    }
}