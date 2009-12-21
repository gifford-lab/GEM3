package edu.mit.csail.cgs.ewok.nouns;

import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

public class CpGIsland extends ScoredRegion {
    private int cpgnum, gcnum, percpg, pergc;
    private float obsexp;

    public CpGIsland(Genome g, String c, int start, int end,
                     int cpgnum, int gcnum, int percpg, int pergc, float obsexp) {
        super(g,c,start,end,0.0);
        this.cpgnum = cpgnum;
        this.gcnum = gcnum;
        this.percpg = percpg;
        this.pergc = pergc;
        this.obsexp = obsexp;
    }
    public int getCpGNum() {return cpgnum;}
    public int getGCNum() {return gcnum;}
    public int getPerCpG() {return percpg;}
    public int getPerGC() {return pergc;}
    public float getObsExp() {return obsexp;}
    public double getScore() {return ((double)percpg) / 100.0;}
}
