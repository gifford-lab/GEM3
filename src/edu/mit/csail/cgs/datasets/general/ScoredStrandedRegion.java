package edu.mit.csail.cgs.datasets.general;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import edu.mit.csail.cgs.datasets.species.Genome;

public class ScoredStrandedRegion extends ScoredRegion implements Stranded {
    private char strand;
    
    public ScoredStrandedRegion(ScoredStrandedRegion copied) {
        super(copied);
        strand = copied.strand;
    }

    public ScoredStrandedRegion(Genome g, String c, int start, int end, double score, char strand) {
        super(g,c,start,end,score);
        this.strand = strand;
    }
    
    public ScoredStrandedRegion(Genome g, DataInputStream dis) throws IOException { 
        super(g,dis);
        strand = dis.readChar();
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        super.save(dos);
        dos.writeChar(strand);
    }

    public char getStrand() {return strand;}

    public boolean equals(Object o) {
        if (o instanceof ScoredStrandedRegion) {
            ScoredStrandedRegion other = (ScoredStrandedRegion)o;
            return super.equals(other) && other.strand == strand;
        } else {
            return false;
        }
    }
    
    public String toString() { 
        String str = getLocationString() + ":" + strand;
        str += " (" + getScore() + ")";
        return str;
    }

}
