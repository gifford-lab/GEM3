package edu.mit.csail.cgs.clustering.vectorcluster;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;

/**
 * A default implementation of VectorClusterElement which implements a 
 * vector of Doubles for each index, with null's indicating missing values.
 * 
 * @author tdanford
 */
public class DefaultVectorClusterElement implements VectorClusterElement {
	
    private Double[] array;
    private Map<String,String> tags;
	
    public DefaultVectorClusterElement(int dim, double min, double max) { 
        double range = max - min;
        Random r = new Random();
        array = new Double[dim];
        tags = new HashMap<String,String>();
        for(int i = 0; i < array.length; i++) { 
            array[i] = r.nextDouble() * range + min;
        }
    }
	
    public DefaultVectorClusterElement(Double[] a) { 
        array = (Double[])a.clone();
        tags = new HashMap<String,String>();
    }
	
    public DefaultVectorClusterElement(double[] a) { 
        array = new Double[a.length];
        for(int i = 0; i < array.length; i++) { array[i] = a[i]; }
        tags = new HashMap<String,String>();
    }
	
    public DefaultVectorClusterElement(Vector<Double> vals) { 
        array = new Double[vals.size()];
        for(int i =0; i < vals.size(); i++) { array[i] = vals.get(i); }
        tags = new HashMap<String,String>();
    }
	
    public void addTag(String k, String v) { 
        tags.put(k, v);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#dimension()
     */
    public int dimension() {
        return array.length;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#getTag(java.lang.String)
     */
    public String getTag(String k) {
        return tags.get(k);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#getValue(int)
     */
    public double getValue(int i) {
        return array[i];
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#hasTag(java.lang.String)
     */
    public boolean hasTag(String k) {
        return tags.containsKey(k);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#isMissingValue(int)
     */
    public boolean isMissingValue(int i) {
        return array[i] == null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.ClusterElement#numMissingValues()
     */
    public int numMissingValues() {
        int c = 0;
        for(int i = 0; i < array.length; i++) { 
            if(isMissingValue(i)) { c += 1; }
        }
        return c;
    }
	
    private static double eps = 1.0e-99;
	
    private static boolean equal(Double d1, Double d2) { 
        if(d1 == null && d2 == null) { return true; }
        if(d1 == null || d2 == null) { return false; }
        if(Math.abs(d1 - d2) <= eps) { return true; }
        return false;
    }
	
    public int hashCode() { 
        int code = 17;
        for(int i = 0; i < array.length; i++) { 
            if(array[i] != null) { 
                long bits = Double.doubleToLongBits(array[i]);
                code += (int)(bits >> 32); code *= 37;
            }
        }
        return code;
    }
	
    public boolean equals(Object o) { 
        if(!(o instanceof DefaultVectorClusterElement)) { return false; }
        DefaultVectorClusterElement e = (DefaultVectorClusterElement)o;
        if(array.length != e.array.length) { return false; }
        for(int i = 0; i < array.length; i++) { 
            if(!equal(array[i], e.array[i])) { return false; }
        }
        return true;
    }
}
