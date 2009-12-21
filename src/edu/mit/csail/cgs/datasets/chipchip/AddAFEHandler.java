package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.utils.parsing.textfiles.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;

public class AddAFEHandler extends AFEHandler {

    private java.sql.Connection cxn;
    private PreparedStatement insert;
    private int idcol, blockcol, colcol, rowcol, ipcol, wcecol, lrcol;
    private boolean doBackgroundSubtraction;

    public AddAFEHandler(Connection c, PreparedStatement i) {
        super();
        cxn = c;
        insert = i;
        idcol = -1; blockcol = -1; colcol = -1;
        rowcol = -1; ipcol = -1; wcecol = -1; lrcol = -1;
        doBackgroundSubtraction = true;
    }

    public void setBackgroundSubtraction(boolean b) {doBackgroundSubtraction = b;}
    public void setIDCol(int i) {idcol = i;}
    public void setBlockCol(int i) {blockcol = i;}
    public void setColCol(int i) {colcol = i;}
    public void setRowCol(int i) {rowcol = i;}
    public void setIPCol(int i) {ipcol = i;}
    public void setWCECol(int i) {wcecol = i;}
    public void setLRCol(int i) {lrcol = i;}
    public String getIDLabel() {return "ProbeName";}
    public String getBlockLabel() {return "Block";}
    public String getColLabel() {return "Col";}
    public String getRowLabel() {return "Row";}
    public String getIPLabel() {return doBackgroundSubtraction ? "rBGSubSignal" : "rMedianSignal";}
    public String getWCELabel() {return doBackgroundSubtraction ? "gBGSubSignal" : "gMedianSignal";}
    public String getLRLabel() {return "LogRatio";}

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
            } else if (headers[i].equals(getLRLabel())) {
                setLRCol(i);
            } 
        }        
        
        /*
        if (idcol == -1) {throw new RuntimeException("Couldn't get ID in " + joinStringArray(headers));}
        //        if (blockcol == -1) {throw new RuntimeException("Couldn't get Block in " + joinStringArray(headers));}
        if (rowcol == -1) {throw new RuntimeException("Couldn't get Row in " + joinStringArray(headers));}
        if (colcol == -1) {throw new RuntimeException("Couldn't get Col in " + joinStringArray(headers));}        
        if (ipcol == -1) {throw new RuntimeException("Couldn't get IP in " + joinStringArray(headers));}        
        if (wcecol == -1) {throw new RuntimeException("Couldn't get WCE in " + joinStringArray(headers));}        
        if (lrcol == -1) {throw new RuntimeException("Couldn't get LogRatio in " + joinStringArray(headers));}
        */        
    }

/* cols are experiment (already set), probeid, block, row, col, ip, wce, mor */
    public void processBodyLine(String fields[]) {        
        try {
            double ip = 0.001 + Double.parseDouble(fields[ipcol]);
            double wce = 0.001 + Double.parseDouble(fields[wcecol]);
            insert.setString(2,fields[idcol]);
            insert.setInt(3,blockcol == -1 ? 1 : Integer.parseInt(fields[blockcol]));
            insert.setInt(4,Integer.parseInt(fields[rowcol]));
            insert.setInt(5,Integer.parseInt(fields[colcol]));
            insert.setDouble(6,ip);
            insert.setDouble(7,wce);
            insert.setDouble(8,Math.pow(10,Double.parseDouble(fields[lrcol])));
            insert.execute();
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(),e);
        } catch (ArrayIndexOutOfBoundsException e) {
            System.err.println("Skipping a line of input because of ");
            e.printStackTrace();
        }



    }


}
