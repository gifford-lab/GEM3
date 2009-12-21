package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;

import java.util.*;
import java.io.*;
import java.sql.*;

public class SQLBayes  extends SQLGeneric implements ChipChipBayes {
    private String dataminsql, maxpostsql, maxstrsql, maxrngpostsql, maxrngstrsql;
    private ArrayList<Double> strength, posterior, strstd, poststd;
    
    public SQLBayes(String analysis, String version, int genomeid) throws NotFoundException {
        super("bayesanalysis",analysis,version,genomeid);
        String exptidstring = " bayesresults.analysis = ? ";
        datasql = "select position, posterior, strength, posteriorstd, strengthstd from bayesresults where " +
            exptidstring + " and chromosome = ? and position >= ? and position <= ? order by position";
        dataminsql = "select position, posterior, strength, posteriorstd, strengthstd from bayesresults where " +
            exptidstring + " and chromosome = ? and position >= ? and position <= ? " + 
            "and posterior >= ? and strength >= ? order by position";           
        maxpostsql = "select max(posterior) from bayesresults where " + exptidstring;
        maxstrsql = "select max(strength) from bayesresults where " + exptidstring;
        maxrngpostsql = "select max(posterior) from bayesresults where " + exptidstring+
            " and chromosome = ? and position >= ? and position <= ?";
        maxrngstrsql = "select max(maxstrength) from bayesresults where "+ exptidstring+
            " and chromosome = ? and position >= ? and position <= ?";
    }
    public void window(String chrom, int start, int stop, double minconf, double minsize) throws NotFoundException {
        if (sameWindow(chrom,start,stop,2)) {return;}
        try {
            PreparedStatement datastmt = chipcxn.prepareStatement(dataminsql);
            datastmt.setFetchSize(400);
            bindExptParams(datastmt);
            datastmt.setInt(sqlparambase+0,getChromID(chrom));
            datastmt.setInt(sqlparambase+1,start);
            datastmt.setInt(sqlparambase+2,stop);
            datastmt.setDouble(sqlparambase+3,minconf);
            datastmt.setDouble(sqlparambase+4,minsize);
            parseWindow(datastmt.executeQuery());
            datastmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    protected void parseWindow(ResultSet rs) throws SQLException {
        ArrayList<Integer> poslist = new ArrayList<Integer>();
        strength = new ArrayList<Double>();
        posterior = new ArrayList<Double>();
        strstd = new ArrayList<Double>();
        poststd = new ArrayList<Double>();
        int i = 0;
        while (rs.next()) {
            poslist.add(rs.getInt(1));
            posterior.add(rs.getDouble(2));
            strength.add(rs.getDouble(3));
            poststd.add(rs.getDouble(4));
            strstd.add(rs.getDouble(5));
            i++;
        }
        count = i;
        int[] positions = new int[poslist.size()];
        for (i = 0; i < poslist.size(); i++) {positions[i] = poslist.get(i);}
        setPositions(positions);
        rs.close();
    }

    public double getPosterior(int i) {return posterior.get(i);}
    public double getStrength(int i) {return strength.get(i);}
    public double getPosteriorStd(int i) {return poststd.get(i);}
    public double getStrengthStd(int i) {return strstd.get(i);}
    public int getReplicates(int i) {return 1;}
    public double getValue(int i, int j) {return getStrength(i);}

    public double getMaxPosterior() {
        try {
            PreparedStatement maxpoststmt = chipcxn.prepareStatement(maxpostsql);
            bindExptParams(maxpoststmt);
            ResultSet rs = maxpoststmt.executeQuery();
            rs.next();
            double result = rs.getDouble(1);
            rs.close();
            maxpoststmt.close();
            return result;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    public double getMaxStrength() {
        try{
            PreparedStatement maxstrstmt = chipcxn.prepareStatement(maxstrsql);
            bindExptParams(maxstrstmt);
            ResultSet rs = maxstrstmt.executeQuery();
            rs.next();
            double result = rs.getDouble(1);
            rs.close();
            maxstrstmt.close();
            return result;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    public double getMaxPosterior(String chrom, int start, int stop) throws NotFoundException {
        try {
            PreparedStatement maxrngpoststmt = chipcxn.prepareStatement(maxrngpostsql);
            bindExptParams(maxrngpoststmt);
            maxrngpoststmt.setInt(sqlparambase+0,getChromID(chrom));
            maxrngpoststmt.setInt(sqlparambase+1,start);
            maxrngpoststmt.setInt(sqlparambase+2,stop);
            ResultSet rs = maxrngpoststmt.executeQuery();
            rs.next();
            double result = rs.getDouble(1);
            rs.close();
            maxrngpoststmt.close();
            return result;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    public double getMaxStrength(String chrom, int start, int stop) throws NotFoundException {
        try {            
            PreparedStatement maxrngstrstmt = chipcxn.prepareStatement(maxrngstrsql);
            bindExptParams(maxrngstrstmt);
            maxrngstrstmt.setInt(sqlparambase+0,getChromID(chrom));
            maxrngstrstmt.setInt(sqlparambase+1,start);
            maxrngstrstmt.setInt(sqlparambase+2,stop);
            ResultSet rs = maxrngstrstmt.executeQuery();
            rs.next();
            double result = rs.getDouble(1);
            rs.close();
            maxrngstrstmt.close();
            return result;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
}
