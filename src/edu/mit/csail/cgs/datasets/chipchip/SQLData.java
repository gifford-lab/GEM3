package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import java.util.*;
import java.io.*;
import java.sql.*;


/**
 * <code>SQLData</code> provides an interface to the data table: channel intensities and
 * ratios from binding experiments on agilent arrays.
 *
 * @author <a href="mailto:arolfe@mit.edu"></a>
 * @version 1.0
 */
public class SQLData extends SQLGeneric implements ChipChipData {
    private String basicdatasql, dataminsql, datamaxcountsql, maxsql,minsql, maxrngsql, minrngsql, exptidstring, probeidsql;
    private Integer[] indexes;
    private double max;
    private double min;
    private int maxCount, alpha;
    private ArrayList<Double> var, ip, wce, ratio;
    private ArrayList<Integer> ids, expts;
    private ArrayList<Character> strands;

    public void close() { 
        super.close();
    }

    /**
     * unlike other SQL* classes, SQLData may need to retrieve data from multiple
     * experiments at once if an experiment has multiple replicates. 
     */
    public SQLData (String exptname, String exptversion, int genomeid, Set<String> replicates) throws NotFoundException {
        super("experiment",exptname,exptversion, genomeid);
        max = -1;
        min = -1;
        ids = new ArrayList<Integer>();
        try {
            if (replicates != null) {
                PreparedStatement stmt = chipcxn.prepareStatement("select id from experiment where name = ? and species = ? and version = ? and replicate = ? and active = 1");
                stmt.setString(1,exptname);
                stmt.setInt(2,speciesid);
                stmt.setString(3,exptversion);
                Iterator<String> reps = replicates.iterator();
                while (reps.hasNext()) {
                    String rep = reps.next();
                    stmt.setString(4,rep);
                    ResultSet rs = stmt.executeQuery();
                    if (rs.next()) {
                        ids.add(rs.getInt(1));
                    } else {
                        throw new NotFoundException("No exptid for " + exptname + "," + exptversion + "," + speciesid +","+rep);
                    }
                    rs.close();
                    stmt.close();
                }
            } else {
                PreparedStatement stmt = chipcxn.prepareStatement("select id from experiment where name = ? and species = ? and version = ?  and active = 1");
                stmt.setString(1,exptname);
                stmt.setInt(2,speciesid);
                stmt.setString(3,exptversion);
                ResultSet rs = stmt.executeQuery();
                while (rs.next()) {
                    ids.add(rs.getInt(1));
                    //System.err.println("GOT ID " + rs.getInt(1));
                }
                rs.close();
                stmt.close();
            }
        } catch (SQLException ex) {
            throw new NotFoundException("Couldn't get experiment id(s) : " + ex,ex);
        }

        if (ids.size() == 0) {
            throw new NotFoundException("No exptid for " + exptname + "," + exptversion + "," + speciesid + "," + replicates);
        }
        init();
    }

    private void init() {
        if (ids.size() == 1) {
            Iterator<Integer> iditer = ids.iterator();
            exptid = iditer.next();
            exptidstring = " data.experiment = ? ";
            sqlparambase = 2;
        } else {
            StringBuilder sb = new StringBuilder("experiment in (");
            for (int i = 0; i < ids.size() - 1; i++) {
                sb.append(ids.get(i));
                sb.append(",");
            }
            sb.append(ids.get(ids.size() -1));
            sb.append(")");
            exptidstring = sb.toString();
            sqlparambase = 1;
        }
        datasql = "select data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, probelocation.strand, data.experiment " + 
        " from data, probelocation where " + exptidstring + " and data.probe = probelocation.id " +
        " and probelocation.chromosome = ? and probelocation.startpos >= ? and probelocation.startpos <= ? " +
        " order by (probelocation.startpos + probelocation.stoppos)/2, data.experiment";
        probeidsql = "select data.probe, probelocation.startpos, probelocation.stoppos " + 
        " from data, probelocation where " + exptidstring + " and data.probe = probelocation.id " +
        " and probelocation.chromosome = ? and probelocation.startpos >= ? and probelocation.startpos <= ? " +
        " order by (probelocation.startpos + probelocation.stoppos)/2, data.experiment";
        basicdatasql = datasql;
        dataminsql = "select data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, probelocation.strand, data.experiment " + 
        " from data, probelocation where " + exptidstring + " and data.ratio >= ? and data.probe = probelocation.id " +
        " and probelocation.chromosome = ? and probelocation.startpos >= ? and probelocation.startpos <= ? " +
        " order by (probelocation.startpos + probelocation.stoppos)/2, data.experiment";
        datamaxcountsql = "select data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, probelocation.strand, data.experiment " + 
        " from data, probelocation where " + exptidstring + " and probelocation.loccount <= ? and data.probe = probelocation.id " +
        " and probelocation.chromosome = ? and probelocation.startpos >= ? and probelocation.startpos <= ? " +
        " order by (probelocation.startpos + probelocation.stoppos)/2, data.experiment";
        maxsql = "select max(ratio) from data where " + exptidstring;
        minsql = "select min(ratio) from data where " + exptidstring;
        maxrngsql = "select count(*), max(data.ratio) from data, probelocation where "+ exptidstring + 
        " and probelocation.chromosome = ? " + 
        " and data.probe = probelocation.id and probelocation.startpos >= ? " + 
        " and probelocation.startpos <= ?";         
        minrngsql = "select count(*), min(data.ratio) from data, probelocation where "+ exptidstring + 
        " and probelocation.chromosome = ? " + 
        " and data.probe = probelocation.id and probelocation.startpos >= ? " + 
        " and probelocation.startpos <= ?";         
        alpha = 0;
        maxCount = -1;
    }

    public void bindExptParams(PreparedStatement ps) throws SQLException {
        if (sqlparambase == 2) {
            //            System.err.println("Binding " + exptid);
            ps.setInt(1,exptid);
        }
    }
    /* returns an SQL where clause for this experiment or set of experiments.
       You must have included the data table and this must be the first
       where clause */
    public String getExptIDString() {
        return exptidstring;
    }
    /* returns the number of replicates in the dataset. 
       This can't just use the replicates provided as input to the constructor
       because some spots are repeated and we need to account for that here.
       In other words, getReplicates() returns the total number of observations
       at a single position. */
    public int getReplicates(int i) {
        if (i == indexes.length - 1) {
            return ip.size() - indexes[i];
        } else if (i < indexes.length - 1) {
            return indexes[i+1] - indexes[i];
        }  else {
            throw new ArrayIndexOutOfBoundsException("i is " + i + " but indexes.length is " + indexes.length + ".  count=" + count);
        }
    }
    /* returns the max ratio of all observations in this data */
    public double getMax() {
        if (max < 0) {
            try {
                PreparedStatement maxstmt = chipcxn.prepareStatement(maxsql);
                ResultSet rs = maxstmt.executeQuery();
                rs.next();
                max = rs.getDouble(1);
                rs.close();
                maxstmt.close();
            } catch (SQLException ex) {
                throw new DatabaseException(ex.toString(),ex);
            }
        }
        return max;
    }
    /* start and stop are chromosomal positions.  Returns the maximum ratio
     between these chromosomal positions, inclusive*/
    public double getMax(String chrom, int start, int stop) throws NotFoundException {
        if (sameWindow(chrom,start,stop,1)) {
            max = (-1*Double.MAX_VALUE);
            for (int i = 0; i < getCount();i++) {
                for (int j = 0; j < getReplicates(i); j++) {
                    if (getRatio(i,j) > max) {max = getRatio(i,j);}
                }
            }
            return max;
        }
        try {
            PreparedStatement maxrngstmt = chipcxn.prepareStatement(maxrngsql);
            bindExptParams(maxrngstmt);
            maxrngstmt.setInt(sqlparambase+0,getChromID(chrom));
            maxrngstmt.setInt(sqlparambase+1,start);
            maxrngstmt.setInt(sqlparambase+2,stop);
            ResultSet rs = maxrngstmt.executeQuery();
            rs.next();
            double maxratio = rs.getDouble(2);
            rs.close();
            maxrngstmt.close();
            return maxratio;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    /* returns the max ratio of all observations in this data */
    public double getMin() {
        if (min < 0) {
            try {
                PreparedStatement minstmt = chipcxn.prepareStatement(minsql);
                ResultSet rs = minstmt.executeQuery();
                rs.next();
                min = rs.getDouble(1);
                rs.close();
                minstmt.close();
            } catch (SQLException ex) {
                throw new DatabaseException(ex.toString(),ex);
            }
        }
        return min;
    }
    /* start and stop are chromosomal positions.  Returns the minimum ratio
    between these chromosomal positions, inclusive*/
    public double getMin(String chrom, int start, int stop) throws NotFoundException {
        if (sameWindow(chrom,start,stop,1)) {
            double min = Double.MAX_VALUE;
            for (int i = 0; i < getCount();i++) {
                for (int j = 0; j < getReplicates(i); j++) {
                    if (getRatio(i,j) <min) {min = getRatio(i,j);}
                }
            }
            return min;
        }
        try {
            PreparedStatement minrngstmt = chipcxn.prepareStatement(minrngsql);
            bindExptParams(minrngstmt);
            minrngstmt.setInt(sqlparambase+0,getChromID(chrom));
            minrngstmt.setInt(sqlparambase+1,start);
            minrngstmt.setInt(sqlparambase+2,stop);
            ResultSet rs = minrngstmt.executeQuery();
            rs.next();
            double maxratio = rs.getDouble(2);
            rs.close();
            minrngstmt.close();
            return maxratio;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
    }
    
    public Vector<int[]> lookupWindowedIdentifiers(String chrom, int start, int stop) throws NotFoundException { 
        try {
            PreparedStatement datastmt = chipcxn.prepareStatement(probeidsql);
            datastmt.setFetchSize(1000);
            bindExptParams(datastmt);
            datastmt.setInt(sqlparambase+0,getChromID(chrom));
            datastmt.setInt(sqlparambase+1,start); 
            datastmt.setInt(sqlparambase+2,stop); 
            
            Vector<int[]> results = new Vector<int[]>();
            
            ResultSet rs = datastmt.executeQuery();
            while(rs.next()) { 
                int[] pair = new int[2];
                
                int probeID = rs.getInt(1);
                int probestart = rs.getInt(2);
                int probeend = rs.getInt(3);
                int pos = (probestart+probeend)/2;
                
                pair[0] = pos; pair[1] = probeID;
                results.add(pair);
            }
            rs.close();
            
            datastmt.close();
            
            return results;
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }  
    }

    public void window(String chrom, int start, int stop, double minratio) throws NotFoundException {
        boolean windowok = sameWindow(chrom,start,stop,2);
        try {
            PreparedStatement datastmt = chipcxn.prepareStatement(dataminsql);
            datastmt.setFetchSize(1000);
            bindExptParams(datastmt);
            datastmt.setInt(sqlparambase+0,getChromID(chrom));
            datastmt.setInt(sqlparambase+1,start); 
            datastmt.setInt(sqlparambase+2,stop); 
            datastmt.setDouble(sqlparambase+3,minratio);
            parseWindow(datastmt.executeQuery());
            datastmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }  
    }

    public void window(String chrom, int start, int stop, int maxcount) throws NotFoundException {
        boolean windowok = sameWindow(chrom,start,stop,3);
        try {
            PreparedStatement datastmt = chipcxn.prepareStatement(datamaxcountsql);
            datastmt.setFetchSize(1000);
            bindExptParams(datastmt);
            datastmt.setDouble(sqlparambase+0,maxcount);
            datastmt.setInt(sqlparambase+1,getChromID(chrom));
            datastmt.setInt(sqlparambase+2,start); 
            datastmt.setInt(sqlparambase+3,stop); 
            parseWindow(datastmt.executeQuery());
            datastmt.close();
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }  
    }

    protected void parseWindow (ResultSet rs) throws SQLException {
        ip = new ArrayList<Double>();
        wce = new ArrayList<Double>();
        ratio = new ArrayList<Double>();
        strands = new ArrayList<Character>();
        ArrayList<Integer> poslist, indexlist;
        poslist = new ArrayList<Integer>();
        indexlist = new ArrayList<Integer>();
        expts = new ArrayList<Integer>();
        int i = 0, lastpos = -1, rep = 0;
        while (rs.next()) {
            int pos = (rs.getInt(4) + rs.getInt(5))/2;
            if (pos != lastpos) {
                poslist.add(pos);
                indexlist.add(i);
            }
            double v = rs.getDouble(1);
            ip.add(v + alpha);
            v = rs.getDouble(2);
            wce.add(v + alpha);
            ratio.add(rs.getDouble(3));
            strands.add(rs.getString(6).charAt(0));
            expts.add(rs.getInt(7));
            lastpos = pos;
            i++;
        }
        count = i;
        if (alpha > 0) {
            for (int j = 0; j < ip.size(); j++) {
                ratio.set(j,(double)ip.get(j) / (double)wce.get(j));
            }
        }
        rs.close();
        var = new ArrayList<Double>(count);
        for (i = 0; i < var.size(); i++) {var.set(i,-1.0);}
        int[] positions = new int[poslist.size()];
        for (i = 0; i < poslist.size(); i++) {positions[i] = poslist.get(i);}
        setPositions(positions);
        indexes = (Integer[])indexlist.toArray(new Integer[0]);
        count = poslist.size();        
    }
    /* get IP value for index i in this window, replicate j */
    public double getIP(int i, int j) {
        return ip.get(indexes[i] + j);
    }
    public double getWCE(int i, int j) {
        return wce.get(indexes[i] + j);
    }    
    public double getRatio(int i, int j)  {
        //        System.err.println("ratio.length is " + ratio.length + " i=" + i + "  j=" + j + "  indexes[i]=" + indexes[i] + " repl[i]=" + getReplicates(i));
        //        System.err.println("   indexes[i+1]=" +  indexes[i+1]);
        return ratio.get(indexes[i] + j);
    }
    public double getValue(int i, int j) {
        return getRatio(i,j);
    }    
    public double getVar(int i, int j) {
        if (var.get(i) < 0) {
            double sum, avg;
            int k, repl = getReplicates(i);
            sum = 0;
            for (k = 0; k < repl; k++) {
                sum += getRatio(i,k);
            }
            avg = sum / repl;
            sum = 0;
            for (k = 0; k < repl; k++) {
                sum += Math.pow(getRatio(i,k) - avg,2);
            }
            sum /= repl;
            var.set(i,sum);
        }
        return var.get(i);
    }
    public char getStrand(int i, int j) {
        return strands.get(indexes[i] + j);
    }
    public int getExptID(int i, int j) {
        return expts.get(indexes[i] + j);
    }
    public int getMaxCount() {return maxCount;}
    public void setMaxCount(int i) {
        maxCount = i;
        if (maxCount == -1) {
            datasql = basicdatasql;
        } else {
            datasql = "select data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, probelocation.strand, data.experiment " + 
            " from data, probelocation where " + exptidstring + " and probelocation.loccount <= " +
            maxCount + " and data.probe = probelocation.id " +
            " and probelocation.chromosome = ? and probelocation.startpos >= ? and probelocation.startpos <= ? " +
            " order by (probelocation.startpos + probelocation.stoppos)/2, data.experiment";
        }


    }
    public void regularize(int alpha) {
        int diff = alpha - this.alpha;
        if (diff == 0) {
            return;
        }
        this.alpha = alpha;
        if (ip == null) {
            return;
        }
        for (int i = 0; i < ip.size(); i++) {
            double v = ip.get(i);
            ip.set(i,v + diff);
            v = wce.get(i);
            wce.set(i,v + diff);
            ratio.set(i,(double)ip.get(i) / (double)wce.get(i));
        }
    }

    /* returns the total number of observations in these experiments */
    public int getTotalCount() throws SQLException {
        String sql = "select count(*) from data where " + exptidstring;

        PreparedStatement countstatement = chipcxn.prepareStatement(sql);
        bindExptParams(countstatement);
        ResultSet rs = countstatement.executeQuery();
        rs.next();
        int result = rs.getInt(1);
        rs.close();
        countstatement.close();
        return result;
    }
    /* returns a 2xn array of all the values in these experiments.  The first
       column are IP values and the second column is WCE 
     */
    public float[][] getRawValues() throws SQLException {
        int count = getTotalCount();        
        String sql = "select channelone, channeltwo from data where " + exptidstring;
        PreparedStatement alldatastatement = chipcxn.prepareStatement(sql);
        bindExptParams(alldatastatement);
        ResultSet rs = alldatastatement.executeQuery();
        float[][] output = new float[2][count];
        int i = 0;
        while (rs.next()) {
            double a = rs.getDouble(1);
            double b = rs.getDouble(2);
            if (i < count) {
                output[0][i] = rs.getFloat(1);
                output[1][i] = rs.getFloat(2);
                i++;
            }
        }
        rs.close();
        alldatastatement.close();
        return output;
    }
    
    /* returns a 2xn array of all the values in these experiments.  The first
    column are IP values and the second column is WCE 
  */
 public float[][] getRawValuesWithIDs() throws SQLException {
     int count = getTotalCount();        
     String sql = "select channelone, channeltwo, data.probe from data where " + exptidstring;
     PreparedStatement alldatastatement = chipcxn.prepareStatement(sql);
     bindExptParams(alldatastatement);
     ResultSet rs = alldatastatement.executeQuery();
     float[][] output = new float[3][count];
     int i = 0;
     while (rs.next()) {
         double a = rs.getDouble(1);
         double b = rs.getDouble(2);
         double c = rs.getDouble(3);
         if (i < count) {
             output[0][i] = rs.getFloat(1);
             output[1][i] = rs.getFloat(2);
             output[2][i] = rs.getFloat(3);
             i++;
         }
     }
     rs.close();
     alldatastatement.close();
     return output;
 }

}
