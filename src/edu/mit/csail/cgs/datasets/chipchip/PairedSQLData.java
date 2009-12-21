package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import java.util.*;
import java.sql.*;

/* sometimes you want to retrieve data from two experiments matched by probe.  
   For example, you might want to generate a scatterplot showing the Cy5 
   intensity in expt A plotted agains the Cy5 intensity in expt B

   You can use this class when A and B are the same experiment, but that'd
   be less efficient than using SQLData.
*/

public class PairedSQLData implements Closeable {
    /* columns[] and the column constants need to match */
    public static final int CHANNELONE = 0, CHANNELTWO = 1, RATIO = 2;
    public static final String[] columns = {"channel one(Cy5)","channel two(Cy3)","ratio"};
    private static final int hashrange = 1000000000;
    private int exptone, expttwo;
    private PreparedStatement basicstmt, samplestmt, limitstmt;
    private java.sql.Connection corecxn, chipcxn;    

    /* takes the two experiment ids for which to retrieve data.
     */
    public PairedSQLData (int one, int two) throws SQLException {
        exptone = one;
        expttwo = two;
        try {
            corecxn = DatabaseFactory.getConnection("core");
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role core",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to database for role core",ex);
        }
        try {
            chipcxn = DatabaseFactory.getConnection("chipchip");
            basicstmt = chipcxn.prepareStatement("select a.channelone, a.channeltwo, a.ratio, b.channelone, b.channeltwo, b.ratio "
                                                 + " from data a, data b where a.probe = b.probe and a.experiment = ? and " 
                                                 + " b.experiment = ?");
            samplestmt = chipcxn.prepareStatement("select * from (select a.channelone, a.channeltwo, a.ratio, b.channelone, b.channeltwo, b.ratio, dbms_random.value as randval "
                                                 + " from data a, data b where a.probe = b.probe and a.experiment = ? and " 
                                                 + " b.experiment = ? order by dbms_random.value) where randval < ?");
            limitstmt = chipcxn.prepareStatement("select aone, atwo, arat, bone, btwo, brat from (select a.channelone aone, a.channeltwo atwo, a.ratio arat, b.channelone bone, b.channeltwo btwo, b.ratio brat "
                                                 + " from data a, data b where a.probe = b.probe and a.experiment = ? and " 
                                                 + " b.experiment = ? order by dbms_random.value) where rownum < ?");
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to database for role chipchip",ex);
        }
    }

    public void setExpts(int one, int two) {
        exptone = one;
        expttwo = two;
    }

    /* returns a 2xn array with the requested data.  Each of colone and coltwo is one
       of CHANNELONE, CHANNELTWO, or RATIO to specify what data you wnat back from 
       each probe
    */
    public float[][] getData(int colone, int coltwo) throws SQLException {
        basicstmt.setInt(1, exptone);
        basicstmt.setInt(2, expttwo);
        return parse(colone, coltwo, basicstmt.executeQuery());
    }
    /* returns a 2xn array containing a sample of the possible data.  samplefrac determines
       what fraction of the eligible datapoints are returned */
    public float[][] getDataSampled(int colone, int coltwo, double samplefrac) throws SQLException {
        throw new RuntimeException("getDataSampled doesn't actually work because Alex couldn't figure out the right SQL for oracle");
//         samplestmt.setInt(1, exptone);
//         samplestmt.setInt(2, expttwo);
//         samplestmt.setInt(3,(int)System.currentTimeMillis());
//         samplestmt.setInt(4, (int)(hashrange * samplefrac));
//         return parse(colone, coltwo, samplestmt.executeQuery());
    }
    /* returns a 2xlimit array containing a sample of the possible data */
    public float[][] getDataLimited(int colone, int coltwo, int limit) throws SQLException {
        limitstmt.setInt(1, exptone);
        limitstmt.setInt(2, expttwo);
        limitstmt.setInt(3,limit);
        return parse(colone, coltwo, limitstmt.executeQuery());
    }    
    private float[][] parse(int colone, int coltwo, ResultSet rs) throws SQLException {
        ArrayList<Float> valsone = new ArrayList<Float>();
        ArrayList<Float> valstwo = new ArrayList<Float>();
        while (rs.next()) {
            valsone.add(rs.getFloat(colone + 1));  // 1 b/c ResultSet is 1 based
            valstwo.add(rs.getFloat(coltwo + 4)); // 4 b/c ResultSet is 1 based and there were three columns for channel one
        }
        rs.close();
        float[][] output = new float[2][valsone.size()];
        for (int i = 0; i < valsone.size(); i++) {
            output[0][i] = valsone.get(i);
            output[1][i] = valstwo.get(i);
        }
        valsone = null;
        valstwo = null;
        return output;
    }
    public void close() {
        try {
            basicstmt.close();
            DatabaseFactory.freeConnection(corecxn);
            DatabaseFactory.freeConnection(chipcxn);
            corecxn = null;
            chipcxn = null;
            basicstmt = null;
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(), e);
        }
    }
    public boolean isClosed() {
        return (basicstmt == null);
    }
}