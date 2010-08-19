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
        if (length < 0) {
            throw new IllegalArgumentException("length must be positive");
        }

        this.length = length;
    }
    public SingleHit (int chrom, int position, float weight, boolean strand, int length) {
        this.chrom = chrom;
        pos = position;
        this.weight = weight;
        this.strand = strand;
        if (length < 0) {
            throw new IllegalArgumentException("length must be positive");
        }
        if (length < (int)Short.MIN_VALUE || length > (int)Short.MAX_VALUE) {
            throw new IllegalArgumentException("length has to fit into a short.  Quit abusing my laziness");
        }

        this.length = (short)length;
    }
    public String toString () {
        return String.format("chrom %d, pos %d, weight %.2f, %s, len %d",
                             chrom, pos, weight, strand ? "+" : "-", length);
    }
    public boolean equals(Object o) {
        if (o instanceof SingleHit) {
            SingleHit other = (SingleHit)o;
            return chrom == other.chrom && pos == other.pos && (Math.abs(weight - other.weight)<.001) && strand == other.strand;
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
    public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chrom;
		result = prime * result + pos;
		result = prime * result + (strand ? 1231 : 1237);
		return result;
	}
}