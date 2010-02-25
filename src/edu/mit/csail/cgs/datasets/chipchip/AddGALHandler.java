package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.io.parsing.textfiles.*;

/* fields of insert are block, col, row, name, id, type, seq */

public class AddGALHandler extends GALHandler {

    private java.sql.Connection cxn;
    private PreparedStatement insert;
    private int linecount;
    private int idcol, seqcol, blockcol, colcol, rowcol, probenamecol, typecol;
    private Set<String> existingKeys;
    private Map<String,String> probeseqs;

    public AddGALHandler (Connection c, PreparedStatement i, Map<String,String> probeseqs) {
        super();
        cxn = c;
        insert = i;
        linecount = 0;
        idcol = -1; seqcol = -1; blockcol = -1; colcol = -1;
        rowcol = -1; probenamecol = -1; typecol = -1;
        existingKeys = null;
        this.probeseqs = probeseqs;
    }

    public void setExistingKeys(Set<String> s) {existingKeys = s;}
    public void setIDCol(int i) {idcol = i;}
    public void setSeqCol(int i) {seqcol = i;}
    public void setBlockCol(int i) {blockcol = i;}
    public void setColCol(int i) {colcol = i;}
    public void setRowCol(int i) {rowcol = i;}
    public void setProbeNameCol(int i) {probenamecol = i;}
    public void setTypeCol(int i) {typecol = i;}
    public String getIDLabel() {return "ID";}
    public String getSequenceLabel() {return "Sequence";}
    public String getBlockLabel() {return "Block";}
    public String getColLabel() {return "Column";}
    public String getRowLabel() {return "Row";}
    public String getNameLabel() {return "Name";}
    public String getControlTypeLabel() {return "ControlType";}
    public void processLabelLine(String[] headers) {
        super.processLabelLine(headers);
        System.err.println("THERE ARE " + headers.length + " columns in the headers");
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals(getIDLabel())) {setIDCol(i);}
            if (headers[i].equals(getSequenceLabel())) {setSeqCol(i);}
            if (headers[i].equals(getBlockLabel())) {setBlockCol(i);}
            if (headers[i].equals(getColLabel())) {setColCol(i);}
            if (headers[i].equals(getRowLabel())) {setRowCol(i);}
            if (headers[i].equals(getNameLabel())) {setProbeNameCol(i);}
            if (headers[i].equals(getControlTypeLabel())) {setTypeCol(i);}
            System.err.println("\t" + i + "\t" + headers[i]);
        }
        if (idcol == -1) {throw new RuntimeException("Couldn't get ID in " + headers);}
        if (seqcol == -1 && probeseqs == null) {
            throw new RuntimeException("Couldn't get Sequence in " + headers);}
        if (rowcol == -1) {throw new RuntimeException("Couldn't get Row in " + headers);}
        if (colcol == -1) {throw new RuntimeException("Couldn't get Col in " + headers);}        
    }
    public Set<String> queryExistingKeys(int designid, int galfileid) throws SQLException {
        HashSet<String> keys = new HashSet<String>();
        Statement stmt = cxn.createStatement();
        ResultSet rs = stmt.executeQuery("select probeid, blockno, rowno, colno from probedesign where " +
                                         " arraydesign = " + designid + " and galfile = " + galfileid);
        while (rs.next()) {
            keys.add(rs.getString(1) + "__" + rs.getInt(2) + "__" + 
                     rs.getInt(3) + "__" + rs.getInt(4));
        }
        rs.close();
        stmt.close();
        return keys;
    }

    public void processBodyLine(String fields[]) {
        String k = fields[idcol] + "__" + (blockcol == -1 ? "1" : fields[blockcol]) + "__" + fields[rowcol] + "__" + fields[colcol];
        if (existingKeys != null && existingKeys.contains(k)) {return;}
        if (fields[idcol].length() > 197) {throw new RuntimeException("Probe name too long : " + fields[idcol]);}
        try {
            insert.setInt(1,blockcol == -1 ? 1 : Integer.parseInt(fields[blockcol]));
            insert.setInt(2,Integer.parseInt(fields[colcol]));
            insert.setInt(3,Integer.parseInt(fields[rowcol]));
            insert.setString(4,probenamecol == -1 ? "" : fields[probenamecol]);
            insert.setString(5,fields[idcol]);
            insert.setString(6,typecol == -1 ? "" : fields[typecol]);
            insert.setString(7,seqcol == -1 ? probeseqs.get(fields[idcol]) : fields[seqcol]);
            insert.execute();
            linecount = (linecount + 1) % 1000;
            if (linecount == 0) {
                cxn.commit();
            }
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }


}
