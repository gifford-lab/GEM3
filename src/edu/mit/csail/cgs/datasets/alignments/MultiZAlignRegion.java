package edu.mit.csail.cgs.datasets.alignments;

import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class MultiZAlignRegion extends ScoredStrandedRegion {

    private Genome otherGenome;
    private String otherChrom;
    private int otherStart, otherEnd;

    public MultiZAlignRegion(Genome baseGenome, String baseChrom, int baseStart, int baseEnd,
                             Genome otherGenome, String otherChrom, int otherStart, int otherEnd,
                             double score, char strand) {
        super(baseGenome,baseChrom,baseStart,baseEnd,score,strand);
        this.otherGenome = otherGenome;
        this.otherChrom = otherChrom;
        this.otherStart = otherStart;
        this.otherEnd = otherEnd;
    }

    public Genome getOtherGenome () {return otherGenome;}
    public String getOtherChrom() {return otherChrom;}
    public int getOtherStart() {return otherStart;}
    public int getOtherEnd() {return otherEnd;}
    public int getOtherWidth() {return otherEnd - otherStart;}
    public String toString() {
        return getGenome().getVersion() + " " + getChrom() + ":" + getStart() + "-" + getEnd() + 
            " to " +
            otherGenome.getVersion() + " " + otherChrom + ":" + otherStart + "-" + otherEnd;
    }
}
