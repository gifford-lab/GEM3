package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class Point implements Comparable<Point> {
	
    /**
     * The genome that this point corresponds to
     */
    private Genome g;
    
    /**
     * The chromosome that this point corresponds to
     */
    private String chrom;
    
    /**
     * The location that this point corresponds to
     */
    private int location;
    
    
    public Point (Genome g, String c, int position) {
        this.g = g;
        chrom = c;
        this.location = position;
    }
    public Point(Point copied){
    	this.g = copied.getGenome();
    	this.chrom = copied.getChrom();
    	this.location = copied.getLocation();
    }
    public Genome getGenome() {return g;}
    public String getChrom() {return chrom;}
    public int getLocation() {return location;}

    public Point clone() {
        return new Point(g,chrom,location);
    }
    
    /**
     * The <code>expand</code> method returns the region which is around <code>distance</code> genomic 
     * coordinates (in bp) from this point.<br>
     * <i>Note</i>: The expansion happens in both directions.<br>
     * However, if <code>distance</code> exceeds any of chromosome ends, the expansion takes place
     * until this end.<br>
     * <b>In this method, indexing starts from 1</b>
     * @param distance The length of bidirectional expansion 
     * @return The region (as a Region) that corresponds to this expansion
     */
    public Region expand(int distance) {
    	int chromLength = g.getChromLength(chrom);
    	
        int ns = location - distance;
        if (ns < 1) {ns = 1;}
        
        int ne = location + distance;
        if (ne > chromLength ) {ne = chromLength;}
        return new Region(g,chrom,ns,ne);
    }
    
    public String toString() {
        return getLocationString();
    }
    
    public String getLocationString() { return chrom + ":" + location; }
    
    public int hashCode() { 
        int code = 17;
        code += g.hashCode(); code *= 37;
        code += chrom.hashCode(); code *= 37;
        code += location; code *= 37;
        return code;
    }
    
    public boolean equals(Object o) {
        if (o instanceof Point) {
            Point r = (Point)o;
            return g.equals(r.getGenome()) &&
                chrom.equals(r.getChrom()) &&
                (location == r.getLocation());
        } else {
            return false;
        }
    }
    
    /** (non-Javadoc)
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    public int compareTo(Point p) {
        if(!chrom.equals(p.chrom)) { return chrom.compareTo(p.chrom); }
        if(location < p.location) { return -1; }
        if(location > p.location) { return 1; }
        return 0;
    }

    //Distance to another point
    public int distance(Point p) { 
        if(!chrom.equals(p.getChrom())) { throw new IllegalArgumentException(p.getChrom()); }
        
        return Math.abs(location-p.getLocation());
    }
    //Offset relative to another point (offset = this-p)
    public int offset(Point p) { 
        if(!chrom.equals(p.getChrom())) { throw new IllegalArgumentException(p.getChrom()); }
        
        return location-p.getLocation();
    }
}
