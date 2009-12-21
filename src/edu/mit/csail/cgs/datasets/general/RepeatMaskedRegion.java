package edu.mit.csail.cgs.datasets.general;

import java.util.*;

import edu.mit.csail.cgs.datasets.species.Genome;

public class RepeatMaskedRegion extends NamedTypedRegion implements Scored {

    private String family;
    private double score;

    public RepeatMaskedRegion (Genome g, String c, int start, int end, 
                               String name, String type, String family,
                               double score, char strand) {
        super(g,c,start,end,name,type,strand);
        this.family = family;
        this.score = score;
    }

    public String getFamily() {return family;}
    public double getScore() {return score / ((getEnd() - getStart() + 1) * 2);}

    public String toString() {
        return "[" + getName() + "," +
            getType() + "," + getFamily() + "]";
    }



}
