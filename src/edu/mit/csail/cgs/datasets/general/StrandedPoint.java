package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class StrandedPoint extends Point implements Stranded {
	
	private char strand;

	public StrandedPoint(Genome g, String c, int start, char str) {
		super(g, c, start);
		strand = str;
	}

	public char getStrand() {
		return strand;
	}

    public int hashCode() { 
        int code = super.hashCode();
        code += strand; code *= 37;
        return code;
    }

    public boolean equals(Object o) { 
        if(!(o instanceof StrandedPoint)) { return false; }
        StrandedPoint sp = (StrandedPoint)o;
        if(strand != sp.strand) { return false; }
        return super.equals(sp);
    }
}
