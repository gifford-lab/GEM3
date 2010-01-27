package edu.mit.csail.cgs.projects.readdb;

public class SingleHit implements Comparable<SingleHit> {
    public int chrom;
    public int pos;
    public float weight;
    public boolean strand;
    public short length;
    public SingleHit (int chrom, int position, float weight, boolean strand, short length) {
        this.chrom = chrom;
        pos = position;
        this.weight = weight;
        this.strand = strand;
        this.length = length;
    }
    public boolean equals(Object o) {
        if (o instanceof HitWeightPair) {
            SingleHit other = (SingleHit)o;
            return chrom == other.chrom && pos == other.pos && weight == other.weight && strand == other.strand;
        } else {
            return false;
        }
    }
    public int compareTo(SingleHit o) {
        if (chrom == o.chrom) {
            return pos - o.pos;
        } else {
            return chrom - o.chrom;
        }
    }
}