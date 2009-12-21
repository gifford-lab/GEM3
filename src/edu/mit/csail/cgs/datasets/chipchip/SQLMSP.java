package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;

import java.util.*;
import java.sql.*;

public class SQLMSP extends SQLGeneric implements ChipChipMSP{
    private ArrayList<Float> ratio, X, pval, pval3, red, green, medianofratios;

    public SQLMSP (String n, String v, int gid) throws NotFoundException {
        super("rosettaanalysis",n,v,gid);
        String exptidstring = " analysis = ? ";
        datasql = "select position, ratio, X, pval, pval3, red, green, medianofratios from rosettaresults " +
            "where " + exptidstring + " and chromosome = ? and position >= ? and position <= ? order by position";
    }
    
    protected void parseWindow(ResultSet rs) throws SQLException {
        ArrayList<Integer> poslist = new ArrayList<Integer>();
        ratio = new ArrayList<Float>();
        X = new ArrayList<Float>();
        pval = new ArrayList<Float>();
        pval3 = new ArrayList<Float>();
        medianofratios = new ArrayList<Float>();
        int i = 0;
        while (rs.next()) {
            poslist.add(rs.getInt(1));
            ratio.add(rs.getFloat(2));
            X.add(rs.getFloat(3));
            pval.add(rs.getFloat(4));
            pval3.add(rs.getFloat(5));
            medianofratios.add(rs.getFloat(8));
            i++;
        }
        count = i;
        int[] positions = new int[poslist.size()];
        for (i = 0; i < poslist.size(); i++) {positions[i] = poslist.get(i);}
        setPositions(positions);
        rs.close();        
    }
    public int getReplicates(int i) {return 1;}
    public float getRatio(int i) {return ratio.get(i);}
    public float getX(int i) {return X.get(i);}
    public double getValue(int i, int j) {
        return getMedianOfRatios(i);
    }
    public float getPval(int i) {return pval.get(i);}
    public float getPval3(int i) {return pval3.get(i);}
    public float getMedianOfRatios(int i) {return medianofratios.get(i);}
    
    public void close() { 
        super.close();
    }

}
