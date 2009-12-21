package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;

public class WeightMatrixSearch {
    public static String getmatrices = "select id from weightmatrix";   
    public static String getmatrixname = "select name, version, type, species from weightmatrix where id = ?";
    public static void main(String args[]) {
        String species = null, name = null, version = null;
        double cutoff = Double.NaN;
        int chunksize = -1, maxresults = -1;
        boolean showmatrix = false, quotable = false;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                species = args[++i];
            }
            if (args[i].equals("--name") ||
                args[i].equals("--matrix")) {
                name = args[++i];
            }
            if (args[i].equals("--version")) {
                version = args[++i];
            }
            if (args[i].equals("--cutoff")) {
                cutoff = Double.parseDouble(args[++i]);
            } 
            if (args[i].equals("--chunksize")) {
                chunksize = Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--maxresults")) {
                maxresults = Integer.parseInt(args[++i]);               
            }
            if (args[i].equals("--showmatrix")) {
                showmatrix = true;
            } 
            if (args[i].equals("--quotable")) {
                quotable = true;
            }
        }
        if (species == null || name == null || version == null || Double.isNaN(cutoff)) {
            System.err.println("WeightMatrixSearch --species --name --version --cutoff");
            System.exit(1);
        }
        try {
            Organism o = new Organism(species);
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            java.sql.Connection corecxn =DatabaseFactory.getConnection("core");
            PreparedStatement namestmt = cxn.prepareStatement(getmatrixname);
            PreparedStatement speciesstmt = corecxn.prepareStatement("select name from species where id = ?");
            
            int wmid = WeightMatrix.getWeightMatrixID(o.getDBID(),
                                                            name,
                                                            version);
            WeightMatrix matrix = WeightMatrix.getWeightMatrix(wmid);
            if (chunksize == -1) {
                chunksize = matrix.length();
            }
            int start = 0;
            Map<WeightMatrix,Double> totalresults = new HashMap<WeightMatrix,Double>();
            while (true) {
                int end = start + chunksize;
                if (end > matrix.length()) {
                    break;
                }
                WeightMatrix querymatrix;
                if (chunksize == matrix.length()) {
                    querymatrix = matrix;
                } else {
                    querymatrix = matrix.subMatrix(start,chunksize);
                }
                Map<WeightMatrix,Double> results = searchByDistance(querymatrix,cutoff);
                totalresults.putAll(results);
                start += 1;
            }
            WeightMatrix keys[] = new WeightMatrix[totalresults.keySet().size()];
            keys = totalresults.keySet().toArray(keys);
            Arrays.sort(keys,new MapComparator(totalresults));
            if (maxresults == -1) {
                maxresults = keys.length;
            }
            for (int i = 0; (i < keys.length && i <maxresults); i++) {
                WeightMatrix hit = keys[i];
                if (quotable) {
                    speciesstmt.setInt(1,matrix.speciesid);
                    ResultSet srs = speciesstmt.executeQuery();
                    srs.next();
                    System.out.println("--species \"" + srs.getString(1) +
                                       "\" --name \"" + hit.name+ 
                                       "\" --version \"" + hit.version + "\"\t\t" + totalresults.get(hit));
                    srs.close();
                } else {
                    System.out.println(hit.name + ", " + hit.version + ", " + hit.type +
                                       " : " + totalresults.get(hit));
                }
                if (showmatrix) {
                    System.out.print(WeightMatrix.printMatrix(hit));
                }
            }
            namestmt.close();
            DatabaseFactory.freeConnection(cxn);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        
    }

    public static Map<WeightMatrix,Double> searchByDistance(WeightMatrix query, double cutoff) {
        return search(query,cutoff,new WMDistanceComparator(true,8));
    }

    public static Map<WeightMatrix,Double> searchByNormalizedDistance(WeightMatrix query, double cutoff) {
        return search(query,cutoff,new WMNormalizedDistanceComparator());
    }

    public static Map<WeightMatrix,Double> search(WeightMatrix query, double cutoff, WMComparator comp) {
        HashMap<WeightMatrix,Double> results = new HashMap<WeightMatrix,Double>();
        try {
            java.sql.Connection cxn =DatabaseFactory.getConnection("annotations");
            PreparedStatement stmt = cxn.prepareStatement(getmatrices);
            ResultSet rs = stmt.executeQuery();
            while (rs.next()) {
                int wmid = rs.getInt(1);
                WeightMatrix target;
                try {
                    target = WeightMatrix.getWeightMatrix(wmid);
                } catch (NotFoundException ex) {                    
                    throw new DatabaseException (ex.toString(),ex);
                }
                double score = comp.compare(query,target);
                if (score <= cutoff) {
                    results.put(target,score);
                }                    
            }
            rs.close();
            stmt.close();
            DatabaseFactory.freeConnection(cxn);
        } catch (SQLException ex) {
            throw new DatabaseException(ex.toString(),ex);
        } catch (UnknownRoleException ex) {
            throw new DatabaseException(ex.toString(),ex);
        }
        return results;
    }

}

class MapComparator implements Comparator<WeightMatrix> {
    private Map<WeightMatrix,Double> data;
    public MapComparator(Map<WeightMatrix,Double> d) {
        data = d;
    }
    public int compare (WeightMatrix first, WeightMatrix second) {
        return Double.compare(data.get(first),data.get(second));
    }
    public boolean equals(Object o) {
        return o == this;
    }
}
