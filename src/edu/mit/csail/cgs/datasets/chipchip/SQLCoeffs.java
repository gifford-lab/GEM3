package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

import java.util.*;
import java.io.*;
import java.sql.*;

public class SQLCoeffs implements ChipChipCoeffs {
    double[] coeffs;
    int min, count;

    public SQLCoeffs(String exptname, String version, int speciesid) throws NotFoundException {
        try {
            java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
            Statement stmt = cxn.createStatement();
            ResultSet rs = stmt.executeQuery("select fragdist from experiment where name = '" + exptname +
                                             "' and version = '" + version + "' and species = " + speciesid);
            if (!rs.next()) {
                throw new NotFoundException("Couldn't get experiment " + exptname +","+version+","+
                                            speciesid);
            }
            int id = rs.getInt(1);
            rs = stmt.executeQuery("select count(*), min(distance) from fragdistentry where distribution = "+
                                   id);
            if (!rs.next()) {  throw new NotFoundException("No results for " + id);}
            count = rs.getInt(1);
            min = rs.getInt(2);
            coeffs = new double[count];
            rs = stmt.executeQuery("select distance, value from fragdistentry where distribution = "+id);
            while (rs.next()) {
                coeffs[rs.getInt(1) - min] = rs.getDouble(2);
            }
            DatabaseFactory.freeConnection(cxn);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get SQLCoeffs " + ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't get SQLCoeffs " + ex.toString(),ex);
        }
    }
    public double getCoeff(int dist) {
        dist = Math.abs(dist);
        if ((dist < min) ||
            (dist - min >= count)) {
            return 0;
        } else {
            return coeffs[dist - min];
        }
    }
    public int getMaxDist() {
        int a = Math.abs(min);
        int b = count + min;
        if (a > b) {
            return a;
        } else {
            return b;
        }
    }
}
