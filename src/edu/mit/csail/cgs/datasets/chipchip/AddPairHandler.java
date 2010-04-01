package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.io.parsing.textfiles.*;

public class AddPairHandler extends PairHandler {

    private java.sql.Connection cxn;
    private PreparedStatement insert;
    private int idcol, colcol, rowcol, valcol;

    public AddPairHandler(Connection c, PreparedStatement i) {
        super();
        cxn = c;
        insert = i;
        idcol = -1; colcol = -1;
        rowcol = -1; valcol = -1;
    }

    public void setIDCol(int i) {idcol = i;}
    public void setColCol(int i) {colcol = i;}
    public void setRowCol(int i) {rowcol = i;}
    public void setValCol(int i) {valcol = i;}
    public String getIDLabel() {return "PROBE_ID";}
    public String getColLabel() {return "X";}
    public String getRowLabel() {return "Y";}
    public String getValLabel() {return "PM";}

    public void processLabelLine(String[] headersone, String[] headerstwo) {
        super.processLabelLine(headersone,headerstwo);
        for (int i = 0; i < headersone.length && i < headerstwo.length; i++) {
            if (!headersone[i].equals(headerstwo[i])) {
                throw new RuntimeException("Mismatched headers : " + headersone[i] + " vs " + headerstwo[i]);
            }
            if (headersone[i].equals(getIDLabel())) {
                setIDCol(i);
            } else if (headersone[i].equals(getRowLabel())) {
                setRowCol(i);
            } else if (headersone[i].equals(getColLabel())) {
                setColCol(i);
            } else if (headersone[i].equals(getValLabel())) {
                setValCol(i);
            } 
        }        
        if (idcol == -1) {throw new RuntimeException("Couldn't get ID in " + headersone);}
        if (rowcol == -1) {throw new RuntimeException("Couldn't get Row in " + headersone);}
        if (colcol == -1) {throw new RuntimeException("Couldn't get Col in " + headersone);}        
        if (valcol == -1) {throw new RuntimeException("Couldn't get PM in " + headersone);}        
    }

    /* cols are experiment (already set), probeid, block, row, col, ip, wce, mor */
    public void processBodyLine(String[] fieldsone, String[] fieldstwo) {
        if (!fieldsone[idcol].equals(fieldstwo[idcol])) {
            throw new RuntimeException("Probes don't match " + fieldsone[idcol] + " vs " + fieldstwo[idcol]);
        }

        double ip = 0.001 + Double.parseDouble(fieldsone[valcol]);
        double wce = 0.001 + Double.parseDouble(fieldstwo[valcol]);
        try {
            insert.setString(2,fieldsone[idcol]);
            insert.setInt(3,1);
            insert.setInt(4,Integer.parseInt(fieldsone[rowcol]));
            insert.setInt(5,Integer.parseInt(fieldsone[colcol]));
            insert.setDouble(6,ip);
            insert.setDouble(7,wce);
            insert.setDouble(8,ip/wce);
            insert.execute();
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        }
    }


}
