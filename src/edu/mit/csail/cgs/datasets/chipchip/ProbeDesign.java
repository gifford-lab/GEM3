/*
 * Created on Aug 3, 2006
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;

/**
 * @author tdanford
 */
public class ProbeDesign {
    
    public static Collection<ProbeDesign> loadDBProbeDesigns() throws SQLException { 
        
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        LinkedList<ProbeDesign> pds = new LinkedList<ProbeDesign>();
        java.sql.Connection c = 
            DatabaseFactory.getConnection(ExptLocator.dbRole);
        Statement s = c.createStatement();
        
        ResultSet rs = s.executeQuery("select arraydesign, blockno, colno, rowno, galfile, " +
        "probename, probeid, type, sequence from probedesign");
        
        while(rs.next()) { 
            ProbeDesign pd = new ProbeDesign(rs, loader);
            pds.addLast(pd);
        }
        
        rs.close();
        s.close();
        DatabaseFactory.freeConnection(c);
        return pds;
    }
    
    public static PreparedStatement prepareLoadByID(java.sql.Connection c) 
        throws SQLException { 
        return c.prepareStatement("select arraydesign, blockno, colno, rowno, galfile, " +
                "probename, probeid, type, sequence from probedesign where id=?");
    }
    
    private int dbid;
    private ArrayDesign arraydesign;
    private int blockno, colno, rowno;
    private GALFile galfile; 
    private String probename, probeid, type, sequence;
    
    public ProbeDesign(ResultSet rs, ChipChipMetadataLoader loader) throws SQLException { 
        int dbid = rs.getInt(1);
        int arrayID = rs.getInt(2);
        blockno = rs.getInt(3);
        colno = rs.getInt(4);
        rowno = rs.getInt(5);
        int galID = rs.getInt(6);
        probename = rs.getString(7);
        probeid = rs.getString(8);
        type = rs.getString(9);
        sequence = rs.getString(10);

        arraydesign = loader.loadArrayDesign(arrayID);
        galfile = loader.loadGALFile(galID);
    }
    
    public int getDBID() { return dbid; }
    public ArrayDesign getArrayDesign() { return arraydesign; }
    public GALFile getGALFile() { return galfile; }
    public int getBlockno() { return blockno; }
    public int getRowno() { return rowno; }
    public int getColno() { return colno; }
    public String getProbename() { return probename; }
    public String getProbeid() { return probeid; }
    public String getSequence() { return sequence; }
    public String getType() { return type; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ProbeDesign)) { return false; }
        ProbeDesign pd = (ProbeDesign)o;
        if(dbid != pd.dbid) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = 17;
        code += dbid; code *= 37;
        return code;
    }
    
}
