package edu.mit.csail.cgs.ewok.nouns;

import java.util.*;
import java.util.regex.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class HarbisonRegCodeProbe extends NamedRegion {
    
    public static final int STRONG = 1;
    public static final int WEAK = 0;
    public static String[] strengthStrings = { "0.005", "0.001" };
    private static Pattern factorPattern = Pattern.compile("([\\w\\d]+)_([\\w\\d-]+)");
    
    public static int getStrength(String str) { 
        for(int i = 0; i < strengthStrings.length; i++) { 
            if(str.equals(strengthStrings[i])) { return i; }
        }
        return -1;
    }
    
    private Map<String,Map<String,Integer>> factorCondScore;
    private int size;

    public HarbisonRegCodeProbe(Genome g, String chrom, int start, int end, String name, 
            int tfCount, String tfList, String bindVals) { 
        super(g,chrom,start,end,name);
        factorCondScore = new HashMap<String,Map<String,Integer>>();
        size = tfCount;
        
        String[] tfArray = tfList.split(",");
        String[] bindArray = bindVals.split(",");
        
        if(tfArray.length < tfCount || bindArray.length < tfCount) { 
            throw new IllegalArgumentException(bindVals + " -- " + tfList); 
        }
        
        for(int i = 0; i < tfCount; i++) { 
            int str = getStrength(bindArray[i]);
            Matcher m = factorPattern.matcher(tfArray[i]);
            if(!m.matches()) { throw new IllegalArgumentException(tfArray[i]); }
            String factor = m.group(1), cond = m.group(2);
            
            if(!factorCondScore.containsKey(factor)) { 
                factorCondScore.put(factor, new HashMap<String,Integer>());
            }
            
            if(factorCondScore.get(factor).containsKey(cond)) { 
                throw new IllegalArgumentException(tfArray[i]);
            }
            
            factorCondScore.get(factor).put(cond, str);
        }
    }
    
    public int size() { return size; }
    
    public boolean isBound(String f, String c, int str) { 
        return factorCondScore.get(f).get(c) == str; 
    }
    
    public int getBindingStrength(String f, String c) { 
        return factorCondScore.get(f).get(c);
    }
    
    public boolean hasBinding(String f, String c) { 
        return factorCondScore.containsKey(f) && factorCondScore.get(f).containsKey(c); 
    }
    
    public Set<String> getFactors() { return factorCondScore.keySet(); }
    public Set<String> getConditions(String f) { return factorCondScore.get(f).keySet(); }
    
    public Set<String> getConditions() { 
        HashSet<String> conds = new HashSet<String>();
        for(String f : factorCondScore.keySet()) { 
            conds.addAll(factorCondScore.get(f).keySet());
        }
        return conds;
    }
    
}
