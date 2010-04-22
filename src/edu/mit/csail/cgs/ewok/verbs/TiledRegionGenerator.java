package edu.mit.csail.cgs.ewok.verbs;

import java.sql.*;
import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;



/** maps a Region to a set of sub-regions that are tiled
 * in a particular array design 
 */

public class TiledRegionGenerator<X extends Region> implements Expander<X,Region> {
    private String arrayDesign;
    private int designID;
    private int spacing, mincount, minlen;

    /**
     * @param design name of the array design 
     * @param spacing the maximum spacing in bp between adjacent
     * probes for them to be counted as part of the same tiled region
     * @param mincount the minimum number of probes in a region (with at most <code>spacing</code> bp between adjacent probes)
     * for a tiled region to be included in the output
     */
    public TiledRegionGenerator(String design, int spacing, int mincount) throws NotFoundException { this(design,spacing,mincount,0);}
    public TiledRegionGenerator(String design, int spacing, int mincount, int minlen) throws NotFoundException {
        arrayDesign = design;
        this.spacing = spacing;
        this.mincount = mincount;
        this.minlen = minlen;
        try {
            java.sql.Connection cxn =
                DatabaseFactory.getConnection("chipchip");
            PreparedStatement ps = cxn.prepareStatement("select id from arraydesign where name = ?");
            ps.setString(1,design);
            ResultSet rs = ps.executeQuery();
            if (rs.next()) {
                designID = rs.getInt(1);
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
            } else {
                rs.close();
                ps.close();
                DatabaseFactory.freeConnection(cxn);
                throw new NotFoundException("Couldn't find array design " + design);
            }
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }

    public Iterator<Region> execute(X r) {
        try {
            Genome g = r.getGenome();
            int chromid = g.getChromID(r.getChrom());
            java.sql.Connection cxn =
                DatabaseFactory.getConnection("chipchip");
            PreparedStatement ps = cxn.prepareStatement("select pl.startpos, pl.stoppos from probelocation pl, probedesign pd where" +
                                                        " pl.id = pd.id and pl.chromosome = ? and pl.stoppos >= ? and pl.startpos <= ? and pd.arraydesign = ? " +
                                                        " and pl.loccount = 1 order by pl.startpos");
//             System.err.println("Looking for probes from " + arrayDesign + ", " + designID + " in " + r + ", " + chromid);
//             System.err.println("spacing=" + spacing + " and mincount="+ mincount);
            ps.setInt(1,chromid);
            ps.setInt(2,r.getStart());
            ps.setInt(3,r.getEnd());
            ps.setInt(4,designID);
            ResultSet rs = ps.executeQuery();
            ArrayList<Region> results = new ArrayList<Region>();
            int laststart = -20 * spacing, laststop = -10 * spacing, count = 0;
            boolean clean = false;
            while (rs.next()) {
                int start = Math.min(rs.getInt(1),rs.getInt(2));
                int stop = Math.max(rs.getInt(1),rs.getInt(2));
                count++;
                if (start - laststop > spacing) {
                    if (laststart < 0) {
                        clean = false;
                    } else {
                        if (count >= mincount) {
                        		if((laststop-laststart)>minlen){
                        			Region found = new Region(g,r.getChrom(),laststart,laststop);
                        			results.add(found);
                        		}
                        }
                        clean = true;
                    }
                    count = 0;
                    laststart = start;
                } else {
                    clean = false;
                }
                laststop = stop;
            }
            if (!clean && (laststart > 0)) {
                results.add(new Region(g,r.getChrom(),laststart,laststop));       
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }

}
