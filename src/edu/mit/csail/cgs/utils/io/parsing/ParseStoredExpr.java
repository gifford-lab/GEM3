/*
 * Created on Jun 1, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.util.*;

/**
 * @author tdanford
 */
public class ParseStoredExpr implements ParseExpr {
    
    private String[] orfs;
    private String[] exptNames;
    private Map<String,Map<String,Double>> storedExpr;
    
    public ParseStoredExpr(Map<String,Map<String,Double>> stored) { 
        storedExpr = stored;
        SortedSet<String> exptNameSet = new TreeSet<String>(storedExpr.keySet());
        SortedSet<String> orfNameSet = new TreeSet<String>();
        exptNames = new String[exptNameSet.size()];
        int i = 0;
        for(String exptName : exptNameSet) { 
            exptNames[i++] = exptName;
			System.out.println(i + ": " + exptNames[i-1] + " (" + exptName + ")");
            orfNameSet.addAll(storedExpr.get(exptName).keySet());
        }
        orfs = new String[orfNameSet.size()];
        i = 0;
        for(String orfName : orfNameSet) { orfs[i++] = orfName; }
    }

    public String getExptName(int i) { 
        return exptNames[i];
    }
    
    public Map<String,Double> getWholeExpt(int i) { 
        return storedExpr.get(exptNames[i]);
    }
    
    public boolean exptHasORF(int i, String orf) { 
        return storedExpr.get(exptNames[i]).containsKey(orf);
    }
    
    public double getValue(int i, String orf) { 
        return storedExpr.get(exptNames[i]).get(orf);
    }
    
    public Double[] getValues(int start, int stop, String orf) { 
        Double[] array = new Double[stop-start+1];
        for(int i = 0; i < start-stop+1; i++) { 
            if(exptHasORF(start+i, orf)) { 
                array[i] = new Double(getValue(start+i, orf));
            } else { 
                array[i] = null;
            }
        }
        return array;
    }
    
    public int numExpts() { 
        return exptNames.length;
    }
    
    public int numORFs(int i) { 
        return storedExpr.get(exptNames[i]).keySet().size();
    }
    
    public double getMax(int i) { 
        return findMax(storedExpr.get(exptNames[i]));
    }
    
    public double getMin(int i) { 
        return findMin(storedExpr.get(exptNames[i]));
    }
    
    private double findMax(Map<String,Double> map) { 
        double e = -Double.MAX_VALUE;
        for(String k : map.keySet()) { 
            if(map.get(k) > e) { e = map.get(k); }
        }
        return e;
    }
    
    private double findMin(Map<String,Double> map) { 
        double e = Double.MAX_VALUE;
        for(String k : map.keySet()) { 
            if(map.get(k) < e) { e = map.get(k); }
        }
        return e;        
    }
}
