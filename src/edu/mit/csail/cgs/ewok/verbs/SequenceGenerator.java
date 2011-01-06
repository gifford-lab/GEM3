package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/** 
 * <code>SequenceGenerator</code> maps a Region to the genomic
 * sequence included in that Region.
 */
public class SequenceGenerator<X extends Region> implements Mapper<X,String>, SelfDescribingVerb {

    private Map<Integer,String> cache;

    // no longer used, but kept for compatibility 
    public SequenceGenerator (Genome g) {        
    }
    public SequenceGenerator() {}
    public void useCache(boolean b) {
        if (b && cache == null) {
            cache = new HashMap<Integer,String>();
        } else {
            cache = null;
        }
    }

    public String execute(X region) {
        String result = null;
        String chromname = region.getChrom();
        
        try {
            Genome genome = region.getGenome();
            int chromid = genome.getChromID(chromname);
            if (cache != null) {
                if (!cache.containsKey(chromid)) {
                    java.sql.Connection cxn =
                DatabaseFactory.getConnection("core");
                    PreparedStatement ps = cxn.prepareStatement("select sequence from chromsequence where id = ?");
                    ps.setInt(1,chromid);
                    ResultSet rs = ps.executeQuery();
                    if (rs.next()) {
                        cache.put(chromid,rs.getString(1));
                    }   
                    rs.close();
                    ps.close();
                    DatabaseFactory.freeConnection(cxn);
                }
                /* oops, couldn't load this sequence into cache */
                if (!cache.containsKey(chromid)) {
                    return null;
                }
                result = cache.get(chromid).substring(region.getStart(), region.getEnd() + 1);
            } else {
                java.sql.Connection cxn =
                DatabaseFactory.getConnection("core");
                PreparedStatement ps;
                int start = Math.max(region.getStart() + 1,0);
                if (DatabaseFactory.isOracle(cxn)) {
                    ps = cxn.prepareStatement("select dbms_lob.substr(sequence,?,?) from chromsequence where id = ?");
                    // length comes first, then offset when using dbms_lob.  This is the opposite of the regular SQL substr
                    ps.setInt(2,start);
                    ps.setInt(1,region.getEnd() - region.getStart() + 1);
                    ps.setInt(3,chromid);                   
                } else {
                    ps = cxn.prepareStatement("select substr(sequence,?,?) from chromsequence where id = ?");
                    ps.setInt(1,start);
                    ps.setInt(2,region.getEnd() - region.getStart() + 1);
                    ps.setInt(3,chromid);                   
                }
                ResultSet rs = ps.executeQuery();
                if (rs.next()) {
                    result = rs.getString(1);
                } 
                rs.close();
                ps.close();
                cxn.commit();
                DatabaseFactory.freeConnection(cxn);
            }
        } catch (SQLException ex) {
            ex.printStackTrace();           
        } catch (UnknownRoleException ex) {
            ex.printStackTrace();
            throw new DatabaseException("Couldn't connect to core",ex);
        }
        if (result.length() != region.getWidth()) {
            System.err.println("Wanted " + region + "(" + 
            		region.getWidth() + ") but only got " + result.length());
        }

        return result;
    }

    private static final String[] inputNames = { "Regions" };
    private static final EchoType[] inputTypes = { new ClassType(Region.class) };
    private static final EchoType outputType = new ClassType(String.class);

    public EchoType[] getInputClasses() { return inputTypes; }

    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }    
}
