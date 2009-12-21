/*
 * Created on Apr 6, 2006
 */
package edu.mit.csail.cgs.ewok.nouns;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class BindingSite extends Region {
    
    private int strand;
    private String factor;
    private double score;
    private String word;

    /**
     * @param g
     * @param c
     * @param start
     * @param end
     */
    public BindingSite(Genome g, String c, int start, int end, int str, String f, double sc) {
        super(g, c, start, end);
        strand = str;
        factor = f;
        score = sc;
        word = null;
    }

    /**
     * @param g
     * @param c
     * @param start
     * @param end
     */
    public BindingSite(Genome g, String c, int start, int end, int str, String f, double sc, String w) {
        super(g, c, start, end);
        strand = str;
        factor = f;
        score = sc;
        word = w;
    }

    public String getFactor() { return factor; }
    public int getStrand() { return strand; }
    public double getScore() { return score; }
    public String getWord() { return word; }
    
    public String toString() { 
        String strstr = strand == 1 ? "+" : "-";
        if(strand == 0) { strstr = "?"; }
        return factor + ": " + strstr;
    }
    
    public int hashCode() { 
        int code = super.hashCode();
        code += factor.hashCode(); code *= 37;
        code += strand; code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof BindingSite)) { return false; }
        BindingSite bs = (BindingSite)o;
        if(!super.equals(bs)) { return false; }
        if(!factor.equals(bs.factor)) { return false; }
        if(strand != bs.strand) { return false; }
        return true;
    }

}
