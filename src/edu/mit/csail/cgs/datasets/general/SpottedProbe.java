/*
 * Created on Mar 13, 2008
 */
package edu.mit.csail.cgs.datasets.general;

import java.util.*;
import java.util.regex.*;

import java.sql.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;


public class SpottedProbe extends NamedRegion {
    
    public static Pattern ratioPattern = Pattern.compile("([^=]+)=([^,]+),([^,]+)");
    
    private Map<String,Double> ratios;
    private Map<String,Double> sds;

    public SpottedProbe(Genome g, ResultSet r) throws SQLException { 
        super(g, r.getString(1), r.getInt(2), r.getInt(3), r.getString(4));
        ratios = new TreeMap<String,Double>();
        sds = new TreeMap<String,Double>();
        String ratstr = r.getString(5);
        String[] a = ratstr.split(":");
        
        for(int i = 0; i < a.length; i++) { 
            Matcher m = ratioPattern.matcher(a[i]);
            if(m.matches()) { 
                String f = m.group(1);
                double ratio = Double.parseDouble(m.group(2));
                double sd = Double.parseDouble(m.group(3));
                ratios.put(f, ratio);
                sds.put(f, sd);
            } else { 
                System.err.println(String.format("Invalid Binding String: %s", a[i]));
            }
        }
    }
    
    public boolean containsRatio(String f) { return ratios.containsKey(f); }
    public double getRatio(String f) { return ratios.get(f); }
    public double getSD(String f) { return sds.get(f); }
    public Collection<String> getFactors() { return ratios.keySet(); }
    
    public Set<String> getBoundFactors(double sdThresh) { 
        TreeSet<String> bound = new TreeSet<String>();
        for(String f : sds.keySet()) { 
            if(sds.get(f) >= sdThresh) { 
                bound.add(f);
            }
        }
        return bound;
    }
}
