package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrixScan;

public class MotifScanResultsGenerator implements Expander <Region,ScoredStrandedRegion> {

    private WeightMatrixScan scan;
    private int matrixid, scanid;
    public MotifScanResultsGenerator (WeightMatrixScan s) throws NotFoundException {
        scan = s;
        if (s.hasscandbid) {
            scanid = s.scandbid;
        } else {            
            try {
                java.sql.Connection cxn =
                    DatabaseFactory.getConnection("annotations");
                PreparedStatement ps;                
                ps = cxn.prepareStatement("select id from weightmatrix where name = ? and version = ?");
                ps.setString(1,scan.matrix.name);
                ps.setString(2,scan.matrix.version);
                ResultSet rs = ps.executeQuery();
                if (rs.next()) {
                    matrixid = rs.getInt(1);
                } else {
                    throw new NotFoundException("Can't find " + scan);
                }
                rs.close();
                ps.close();
                ps = cxn.prepareStatement("select id from weightmatrixscan where weightmatrix = ? and name = ?");
                ps.setInt(1,matrixid);
                ps.setString(2,scan.scanname);
                rs = ps.executeQuery();
                if (rs.next()) {
                    scanid = rs.getInt(1);
                } else {
                    throw new NotFoundException("Can't find " + scan);
                }
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
            } catch (SQLException ex) {
                throw new DatabaseException("Couldn't get set up for " + scan,ex);
            } catch (UnknownRoleException ex) {
                throw new DatabaseException("Couldn't connect for role annotations",ex);
            }
        }    
    }

    public Iterator<ScoredStrandedRegion> execute(Iterator<Region> r) {
        return new ExpanderIterator<Region,ScoredStrandedRegion>(this,r);
    }

    public Iterator<ScoredStrandedRegion>execute(Region r) {
        try {
            java.sql.Connection cxn =
                DatabaseFactory.getConnection("annotations");
            String sql = "select startpos,stoppos,strand,score from wms_hits where scan = ? and chromosome = ? and startpos >= ? and stoppos <= ? order by startpos";
            PreparedStatement ps = cxn.prepareStatement(sql);
            ps.setInt(1,scanid);
            ps.setInt(2,r.getGenome().getChromID(r.getChrom()));
            ps.setInt(3,r.getStart());
            ps.setInt(4,r.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<ScoredStrandedRegion> results = new ArrayList<ScoredStrandedRegion>();
            while (rs.next()) {
                results.add(new ScoredStrandedRegion(r.getGenome(),
                                                     r.getChrom(),
                                                     rs.getInt(1),
                                                     rs.getInt(2),
                                                     rs.getDouble(4),
                                                     rs.getString(3).charAt(0)));
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get results for WeightMatrixScan " + scan + " in " + r ,ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect for role annotations",ex);
        }
    }

}
