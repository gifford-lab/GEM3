package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class NamedStrandedRegion extends StrandedRegion implements Named {

    private String name;
    
    public NamedStrandedRegion(NamedStrandedRegion r) { 
    	super(r);
    	name = r.name;
    }
    
    public NamedStrandedRegion(Region r, String n, char st) { 
    	super(r, st);
    	name = n;
    }

    public NamedStrandedRegion(Genome g, String chrom, int start, int end, String name, char s) {
        super(g,chrom,start,end, s);
        this.name = name;
    }
    
    public String getName() {return name;}
    public void setName(String n) { name = n; }
    public String toString() {
        return name + " " + super.toString();
    }
}
