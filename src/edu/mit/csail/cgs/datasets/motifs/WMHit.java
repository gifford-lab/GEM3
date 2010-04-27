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
    public String toString() {
        return String.format("%.2f at %d:%d-%d:%s",
                             score,chromid,start,end,strand);
    }
    public int getStart() {return start;}
    public int getEnd() {return end;}
    public int getScanID() {return scanid;}
    public int getChromID() {return chromid;}
    public String getStrand() {return strand;}
    public float getScore() {return score;}
}
