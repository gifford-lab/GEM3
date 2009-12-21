package edu.mit.csail.cgs.ewok.nouns;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Scored;
import edu.mit.csail.cgs.datasets.species.Genome;

public class HarbisonRegCodeRegion extends NamedRegion implements Scored {
    private double score;
    private String chip;
    private int cons;

    public HarbisonRegCodeRegion(Genome g, String chrom, int start, int end, String name,
                                 double score, String chipEvidence, int consSpecies) {
        super(g,chrom,start,end,name);
        this.score = score;
        chip = chipEvidence;
        cons = consSpecies;
    }
    public double getScore() {return score;}
    public String getChipEvidence() {return chip;}
    public int getConsSpecies() {return cons;}

}
