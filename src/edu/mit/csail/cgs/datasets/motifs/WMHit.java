package edu.mit.csail.cgs.datasets.motifs;

public class WMHit {
    public int scanid, chromid, start, end;
    public String strand;
    public float score;
    public WMHit (int scanid, int chromid, int start, int end, String strand, float score) {
        this.scanid = scanid;
        this.chromid = chromid;
        this.start = start;
        this.end = end;
        this.strand = strand;
        this.score = score;
    }
}
