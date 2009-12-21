package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;

import java.util.*;
import java.io.*;
import java.sql.*;

public class SQLMLE extends SQLGeneric implements ChipChipMLE {
    private int lasttype;
    private String dataminsql, maxsql,minsql, maxrngsql, minrngsql;
    private ArrayList<Double> b_i, bindll, nullll, lograt, conf;

    public SQLMLE (String analysis, String version, int genomeid)  
        throws NotFoundException {
        super("mleanalysis",analysis,version,genomeid);
        String exptidstring = " analysis = ? ";
        datasql = "select position, b_i, bindll, nullll, lograt, conf " +
            "from mleresults where " + exptidstring + " and chromosome = ?" +
            " and position >= ? and position <= ? order by position";
        dataminsql = "select position, b_i, bindll, nullll, lograt, conf " +
            "from mleresults where " + exptidstring + " and chromosome = ?" +
            " and position >= ? and position <= ? " +
            " and b_i >= ? and conf <= ? order by position";
        maxsql = "select max(b_i) from mleresults where " + exptidstring;
        minsql = "select min(b_i) from mleresults where " + exptidstring;
        maxrngsql = "select max(b_i) from mleresults where " +
            exptidstring + " and chromosome = ? " + 
            " and position >= ? and position <= ? ";
        minrngsql = "select min(b_i) from mleresults where " +
        	exptidstring + " and chromosome = ? " + 
        	" and position >= ? and position <= ? ";
    }
    public int getReplicates(int i) {
        return 1;
    }
    public double getMax() {
        try {
            PreparedStatement maxstmt = chipcxn.prepareStatement(maxsql);
            maxstmt.setInt(1,exptid);
            ResultSet rs = maxstmt.executeQuery();
            if(rs.next()) {
                double d = rs.getDouble(1);
                rs.close();
                maxstmt.close();
                return d;
            } else {
                throw new DatabaseException("Can't get max for analysis " + exptid);
            }
        }catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    /* start and stop are chromosomal positions */
    public double getMax(String chrom, int start, int stop) throws NotFoundException {
        try {
            PreparedStatement maxrngstmt = chipcxn.prepareStatement(maxrngsql);
            bindExptParams(maxrngstmt);
            maxrngstmt.setInt(sqlparambase+0,exptid);
            maxrngstmt.setInt(sqlparambase+1,getChromID(chrom));
            maxrngstmt.setInt(sqlparambase+2,start);
            maxrngstmt.setInt(sqlparambase+3,stop);
            ResultSet rs = maxrngstmt.executeQuery();
            rs.next();
            double d = rs.getDouble(1);
            rs.close();
            maxrngstmt.close();
            return d;
        }catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    public double getMin() {
        try {
            PreparedStatement minstmt = chipcxn.prepareStatement(minsql);
            minstmt.setInt(1,exptid);
            ResultSet rs = minstmt.executeQuery();
            if(rs.next()) {
                double d = rs.getDouble(1);
                rs.close();
                minstmt.close();
                return d;
            } else {
                throw new DatabaseException("Can't get min for analysis " + exptid);
            }
        }catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    /* start and stop are chromosomal positions */
    public double getMin(String chrom, int start, int stop) throws NotFoundException {
        try {
            PreparedStatement minrngstmt = chipcxn.prepareStatement(minrngsql);
            bindExptParams(minrngstmt);
            minrngstmt.setInt(sqlparambase+0,exptid);
            minrngstmt.setInt(sqlparambase+1,getChromID(chrom));
            minrngstmt.setInt(sqlparambase+2,start);
            minrngstmt.setInt(sqlparambase+3,stop);
            ResultSet rs = minrngstmt.executeQuery();
            rs.next();
            double d = rs.getDouble(1);
            rs.close();
            minrngstmt.close();
            return d;
        }catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    public void window(String chrom, int start, int stop, double minsize, double maxconf) 
    throws NotFoundException {
        if (sameWindow(chrom,start,stop,2)) {return;}
        try {
            PreparedStatement datastmt = chipcxn.prepareStatement(dataminsql);     
            datastmt.setFetchSize(100);
            bindExptParams(datastmt);
            datastmt.setInt(sqlparambase+0,getChromID(chrom));
            datastmt.setInt(sqlparambase+1,start);
            datastmt.setInt(sqlparambase+2,stop);
            datastmt.setDouble(sqlparambase+3,minsize);
            datastmt.setDouble(sqlparambase+4,maxconf);
            parseWindow(datastmt.executeQuery());
            datastmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    
    protected void parseWindow(ResultSet rs) throws SQLException {
        ArrayList<Integer> poslist = new ArrayList<Integer>();
        b_i = new ArrayList<Double>();
        bindll = new ArrayList<Double>();
        nullll = new ArrayList<Double>();
        lograt = new ArrayList<Double>();
        conf = new ArrayList<Double>();
        int i = 0;
        while (rs.next()) {
            poslist.add(rs.getInt(1));
            b_i.add(rs.getDouble(2));
            bindll.add(rs.getDouble(3));
            nullll.add(rs.getDouble(4));
            lograt.add(rs.getDouble(5));
            conf.add(rs.getDouble(6));
            i++;
        }
        count = i;
        int[] positions = new int[poslist.size()];
        for (i = 0; i < poslist.size(); i++) {positions[i] = poslist.get(i);}
        setPositions(positions);
        rs.close();
    }

    public double getSize(int i) {
        return b_i.get(i);
    }
    public double getValue(int i, int j) {
        return getSize(i);
    }
    public double getBindLL(int i){
        return bindll.get(i);
    }
    public double getNullLL(int i){
        return nullll.get(i);
    }
    public double getLogRat(int i){
        return lograt.get(i);
    }
    public double getConf(int i) {
        return conf.get(i);
    }
    public void close() {
        super.close();
    }


}
