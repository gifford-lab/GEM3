/*
 * Created on Sep 6, 2006
 */
package edu.mit.csail.cgs.utils.io.parsing.affyexpr;

/**
 * @author tdanford
 */
public class AffyProbe {
    
    private String probeName;
    private String geneSymbol, unigene, locusID, geneName;

    public AffyProbe(String pn, String gs, String ug, String li, String gn) {
        probeName = pn;
        geneSymbol = gs;
        unigene = ug;
        locusID = li;
        geneName = gn;
    }
    
    public AffyProbe(String line) { 
        String[] array = line.split("\t");
        probeName = array[0];
        geneSymbol = array[1];
        unigene = array[2];
        locusID = array[3];
        geneName = array[4];
    }
    
    public String getProbeName() { return probeName; }
    public String getGeneSymbol() { return geneSymbol; }
    public String getUnigene() { return unigene; }
    public String getLocusID() { return locusID; }
    public String getGeneName() { return geneName; }
    
    public boolean equals(Object o) { 
        if(!(o instanceof AffyProbe)) { return false; }
        AffyProbe ap = (AffyProbe)o;
        return probeName.equals(ap.probeName);
    }
    
    public int hashCode() { return probeName.hashCode(); }
    public String toString() { return probeName; }

}
