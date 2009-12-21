/*
 * Created on Apr 18, 2006
 */
package edu.mit.csail.cgs.ewok.nouns;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Point;

/**
 * @author tdanford
 */
public class ScoredPoint extends Point {
    
    private Double motifScore, bgScore;

    public ScoredPoint(Genome g, String c, int start, Double ms, Double bs) {
        super(g, c, start);
        motifScore = ms; bgScore = bs;
    }
    
    public ScoredPoint(ScoredPoint target) { 
        super(target.getGenome(), target.getChrom(), target.getLocation());
        motifScore = target.motifScore;
        bgScore = target.bgScore;
    }
    
    public double getMotifScore() { return motifScore; }
    public double getBgScore() { return bgScore; }
    public boolean hasMotifScore() { return motifScore != null; }
    public boolean hasBgScore() { return bgScore != null; }
    
    public double getRatio() { return motifScore / bgScore; }
    public double getDifference() { return motifScore - bgScore; }
    
    public ScoredPoint addScoredPoint(ScoredPoint sp) {
        ScoredPoint ns = new ScoredPoint(this);
        ns.motifScore = addValues(ns.motifScore, sp.motifScore);
        ns.bgScore = addValues(ns.bgScore, sp.bgScore);
        return ns;
    }
    
    public ScoredPoint multiplyScoredPoint(ScoredPoint sp) { 
        ScoredPoint ns = new ScoredPoint(this);
        ns.motifScore = multiplyValues(ns.motifScore, sp.motifScore);
        ns.bgScore = multiplyValues(ns.bgScore, sp.bgScore);
        return ns;        
    }
    
    private static Double addValues(Double v1, Double v2) { 
        if(v1 == null && v2 == null) { return null; }
        if(v1 == null) { return v2; }
        if(v2 == null) { return v1; }
        return v1 + v2;
    }

    private static Double multiplyValues(Double v1, Double v2) { 
        if(v1 == null && v2 == null) { return null; }
        if(v1 == null) { return v2; }
        if(v2 == null) { return v1; }
        return v1 * v2;
    }

    /* (non-Javadoc)
     * @see java.lang.Comparable#compareTo(T)
     */
    public int compareTo(ScoredPoint ss) {
        return super.compareTo(ss);
    }
    
    public int hashCode() { 
        int code = super.hashCode();

        long bits = 0;
        
        if(motifScore != null) { 
            bits = Double.doubleToLongBits(motifScore);
            code += (int)(bits >> 32); code *= 37;
        }
        
        if(bgScore != null) { 
            bits = Double.doubleToLongBits(bgScore);
            code += (int)(bits >> 32); code *= 37;
        }
        
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof ScoredPoint)) { return false; }
        ScoredPoint ss = (ScoredPoint)o;
        if(!super.equals(ss)) { return false; }
        
        if(motifScore != null || ss.motifScore != null) { 
            if(motifScore == null || ss.motifScore == null) { return false; }
            if(!motifScore.equals(ss.motifScore)) { return false; }
        }
        
        if(bgScore != null || ss.bgScore != null) { 
            if(bgScore == null || ss.bgScore == null) { return false; }
            if(!bgScore.equals(ss.bgScore)) { return false; }
        }
        
        return true;
    }

}
