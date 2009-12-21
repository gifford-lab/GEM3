/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.io.*;
import java.text.*;
import java.util.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Point;

/**
 * @author tdanford
 */
public class Probe extends Point {
    
    private Map<String,double[]> values;

    public Probe(Genome genome, String chrom, int loc) {
        super(genome, chrom, loc);
        values = new HashMap<String,double[]>();
    }

    public boolean containsKey(String k) { return values.containsKey(k); }
    public void addValue(String k, double[] val) { values.put(k, val); }
    
    public double getAbsoluteMaxValue() { 
        double max = 0.0;
        for(String k : values.keySet()) {
            double[] array = values.get(k);
            for(int i = 0; i < array.length; i++) {
                if(array[i] >= max) { 
                    max = array[i];
                }
            }
        }
        
        return max;
    }
    
    public double getTagMaxValue(String tag) { 
        double max = 0.0;
        for(String k : values.keySet()) {
            if(k.endsWith(tag)) { 
                double[] array = values.get(k);
                for(int i = 0; i < array.length; i++) {
                    if(array[i] >= max) { 
                        max = array[i];
                    }
                }
            }
        }
        
        return max;        
    }
    
    public double getMaxValue(String k) { 
        double max = 0.0;
        double[] array = values.get(k);
        for(int i = 0; i < array.length; i++) { 
            if(i == 0 || array[i] > max) { 
                max = array[i];
            }
        }
        return max;
    }
    
    public void addValue(String k, double v) { 
        double[] vals = new double[1]; vals[0] = v;
        values.put(k, vals);
    }
    
    public double[] getValue(String k) { return values.get(k); }    
    public Collection<String> getKeys() { return values.keySet(); }
    public double getAverageValue(String k) {
        double sum = 0;
        double[] array = values.get(k);
        if (array == null) {
            System.err.println("KEYS are " + values.keySet());
            throw new NullPointerException("NO SUCH KEY " + k);
        }
        for (int i = 0; i < array.length; i++) {
            sum += array[i];
        }
        return sum / array.length;
    }
    
    public void debugPrint(PrintStream ps) {
        NumberFormat nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        
        ps.println("Probe (" + getLocationString() + ")");
        for(String key : values.keySet()) { 
            double[] array = values.get(key);
            ps.print("\t" + key + " : ");
            for(int i = 0; i < array.length; i++) { 
                if(i > 0) { ps.print(", "); }
                ps.print(nf.format(array[i]));
            }
            ps.println();
        }
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append(super.toString());
        for(String key : values.keySet()) { 
            sb.append(" " + key + ":");
			for(int i = 0; i < values.get(key).length; i++) {
				sb.append(" " + values.get(key)[i]);
			}
        }
        return sb.toString();
    }

    public boolean equals(Object o) { 
        if(!(o instanceof Probe)) { return false; }
        return super.equals((Probe)o);
    }
    
    public void debugPrintValueKeys() {
        int i = 1;
        for(String k : values.keySet()) { 
            System.err.println(i + ": " + k);
            i++;
        }
    }
}
