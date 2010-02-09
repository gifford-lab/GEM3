package edu.mit.csail.cgs.datasets.chipchip;

import java.sql.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.io.parsing.textfiles.*;

public class AddGPRHandler extends GPRHandler {

    private java.sql.Connection cxn;
    private PreparedStatement insert;
    private int idcol, blockcol, colcol, rowcol, ipcol, wcecol, morcol;

    public AddGPRHandler(Connection c, PreparedStatement i) {
        super();
        cxn = c;
        insert = i;
        idcol = -1; blockcol = -1; colcol = -1;
        rowcol = -1; ipcol = -1; wcecol = -1; morcol = -1;
    }

    public void setIDCol(int i) {idcol = i;}
    public void setBlockCol(int i) {blockcol = i;}
    public void setColCol(int i) {colcol = i;}
    public void setRowCol(int i) {rowcol = i;}
    public void setIPCol(int i) {ipcol = i;}
    public void setWCECol(int i) {wcecol = i;}
    public void setMORCol(int i) {morcol = i;}
    public String getIDLabel() {return "ID";}
    public String getBlockLabel() {return "Block";}
    public String getColLabel() {return "Column";}
    public String getRowLabel() {return "Row";}
    public String getIPLabel() {return "F635 Median - B635";}
    public String getWCELabel() {return "F532 Median - B532";}
    public String getMORLabel() {return "Median of Ratios";}

    public void processLabelLine(String[] headers) {
        super.processLabelLine(headers);
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals(getIDLabel())) {
                setIDCol(i);
            } else if (headers[i].equals(getBlockLabel())) {
                setBlockCol(i);
            } else if (headers[i].equals(getRowLabel())) {
                setRowCol(i);
            } else if (headers[i].equals(getColLabel())) {
                setColCol(i);
            } else if (headers[i].equals(getIPLabel())) {
                setIPCol(i);
            } else if (headers[i].equals(getWCELabel())) {
                setWCECol(i);
            } else if (headers[i].equals(getMORLabel())) {
                setMORCol(i);
            } 
        }        
        String headerString = "";
        for (int i = 0; i < headers.length; i++) {
            headerString = headerString + "\t" + headers[i];
        }

        if (idcol == -1) {throw new RuntimeException("Couldn't get ID in " + headerString);}
        if (blockcol == -1) {throw new RuntimeException("Couldn't get Block in " + headerString);}
        if (rowcol == -1) {throw new RuntimeException("Couldn't get Row in " + headerString);}
        if (colcol == -1) {throw new RuntimeException("Couldn't get Col in " + headerString);}        
        if (ipcol == -1) {throw new RuntimeException("Couldn't get IP in " + headerString);}        
        if (wcecol == -1) {throw new RuntimeException("Couldn't get WCE in " + headerString);}        
        if (morcol == -1) {throw new RuntimeException("Couldn't get MOR in " + headerString);}        
    }

/* cols are experiment (already set), probeid, block, row, col, ip, wce, mor */
    public void processBodyLine(String fields[]) {
        double ip = 0.001 + Double.parseDouble(fields[ipcol]);
        double wce = 0.001 + Double.parseDouble(fields[wcecol]);
        try {
            insert.setString(2,fields[idcol]);
            insert.setInt(3,Integer.parseInt(fields[blockcol]));
            insert.setInt(4,Integer.parseInt(fields[rowcol]));
            insert.setInt(5,Integer.parseInt(fields[colcol]));
            insert.setDouble(6,ip);
            insert.setDouble(7,wce);
            insert.setDouble(8,Double.parseDouble(fields[morcol]));
            insert.execute();
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        }
    }


}
