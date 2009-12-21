package edu.mit.csail.cgs.datasets.chipchip;

import java.util.Set;
import java.util.HashSet;
import java.util.Map;
import java.sql.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.parsing.textfiles.*;

public class AddTDTHandler extends RowsColumnsHandler {

    private Connection cxn;
    private PreparedStatement insert;
    private int colcol, rowcol, namecol, idcol, typecol, seqcol;
    private Set<String> existingKeys;
    private Map<String,String> probeseqs;

    public AddTDTHandler(Connection cxn,
                         PreparedStatement insert,
                         Map<String,String> probeseqs) {
        this.cxn = cxn;
        this.insert = insert;
        this.probeseqs = probeseqs;
    }

    public void processLabelLine(String[] fields) {
        super.processLabelLine(fields);
        try {
            rowcol = getColumn("Row");
        } catch (NullPointerException e) {
            throw new RuntimeException("Can't find Row column");
        }
        try {
            colcol = getColumn("Column");
        } catch (NullPointerException e) {
            throw new RuntimeException("Can't find Column column");
        }
        try {
            idcol = getColumn("ID");
        } catch (NullPointerException e) {
            throw new RuntimeException("Can't find ID column");
        }
        try {
            seqcol = getColumn("Sequence");
        } catch (NullPointerException e) {
            if (probeseqs == null) {
                throw new RuntimeException("Can't find Sequence column");
            } else {
                seqcol = -1;
            }
        }
        try {
            namecol = getColumn("Name");
        } catch (NullPointerException e) {
            namecol = -1;
        }
        try {
            typecol = getColumn("ControlType");
        } catch (NullPointerException e) {
            typecol = -1;
        }

    }

    /* insert needs
       block, col, row, probename, probeid, type, sequence */
    public void processBodyLine(String[] fields) {
        try {
            String k = fields[idcol] + "__1__" + (fields[rowcol] + "__" + fields[colcol]);
            if (existingKeys != null && existingKeys.contains(k)) {return;}
            insert.setInt(1,1);
            insert.setInt(2,Integer.parseInt(fields[colcol]));
            insert.setInt(3,Integer.parseInt(fields[rowcol]));
            if (namecol != -1) {
                insert.setString(4,fields[namecol]);
            } else {
                insert.setString(4,null);
            }
            insert.setString(5,fields[idcol]);
            if (typecol != -1) {
                insert.setString(6,fields[typecol]);
            } else {
                insert.setString(6,"");
            }
            insert.setString(7,seqcol == -1 ? probeseqs.get(fields[idcol]) : fields[seqcol]);
            insert.execute();

        } catch (SQLException e) {
            e.printStackTrace();
            throw new DatabaseException(e.toString(),e);
        } catch (ArrayIndexOutOfBoundsException e) {
            e.printStackTrace();
            System.err.println("fields.length is " + fields.length);
            for (int i = 0; i < fields.length; i++) {
                System.err.println("  i=" + i +"  " + fields[i]);
            }

        } catch (NumberFormatException e) {
            e.printStackTrace();
            for (int i = 0; i < fields.length; i++) {
                System.err.println("  i=" + i +"  " + fields[i]);
            }

        } 

    }

    public void setExistingKeys(Set<String> s) {existingKeys = s;}

    public Set<String> queryExistingKeys(int designid, int galfileid) throws SQLException {
        HashSet<String> keys = new HashSet<String>();
        Statement stmt = cxn.createStatement();
        ResultSet rs = stmt.executeQuery("select probeid, blockno, rowno, colno from probedesign where " +
                                         " arraydesign = " + designid + " and galfile = " + galfileid);
        while (rs.next()) {
            keys.add(rs.getString(1) + "__" + rs.getInt(2) + "__" + 
                     rs.getInt(3) + "__" + rs.getInt(4));
        }
        System.err.println("Found " + keys.size() + " existing keys");
        rs.close();
        stmt.close();
        return keys;
    }
}

