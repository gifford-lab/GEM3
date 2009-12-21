/*
 * Created on Jul 16, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.hypotheses.database;

import java.sql.*;
import java.util.*;

import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.*;
import java.sql.Connection;

public class HypothesisLoader implements Closeable {
    
    public static void main(String[] args) { 
        try {
            HypothesisLoader loader = new HypothesisLoader();
            
            for(Investigation inv : loader.loadInvestigations()) { 
                System.out.println(inv);
            }
            
            loader.close();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }
    
    public static final String role = "hypothesis";

    private Connection cxn;
    
    public HypothesisLoader() throws SQLException { 
        cxn = DatabaseFactory.getConnection(role);
    }
    
    public Set<Investigation> loadInvestigations() throws SQLException {
        HashSet<Investigation> names = new HashSet<Investigation>();
        PreparedStatement ps = cxn.prepareStatement("select id, name from investigation");
        ResultSet rs = ps.executeQuery();
        while(rs.next()) { 
            Investigation invest = new Investigation(rs);
            names.add(invest);
        }
        rs.close();
        ps.close();
        return names;
    }
    
    public Investigation getInvestigation(String invest) throws SQLException { 
        PreparedStatement ps = cxn.prepareStatement("select id, name from investigation where name=?");
        ps.setString(1, invest);
        ResultSet rs = ps.executeQuery();
        Investigation inv = null;
        
        if(rs.next()) { 
            inv = new Investigation(rs);
        }
        rs.close();
        ps.close();

        if(inv == null) { 
            String seq = Sequence.getInsertSQL(cxn, "id");
            ps = cxn.prepareStatement("insert into investigation (id, name) values (" + seq + ", ?)");
            ps.setString(1, invest);
            ps.executeUpdate();
            ps.close();

            return getInvestigation(invest);
        } else { 
            return inv;
        }
    }

    public void close() {
        DatabaseFactory.freeConnection(cxn);
        cxn = null;
    }

    public boolean isClosed() {
        return cxn == null;
    }
}
