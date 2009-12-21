package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class ScoredPoint extends Point implements Scored {
	
	private double score;
	
    public ScoredPoint (Genome g, String c, int position, double s) {
    	super(g, c, position);
    	score = s;
    }
    
    public double getScore() { return score; }
    
    public ScoredPoint clone() {
        return new ScoredPoint(getGenome(), getChrom(), getLocation(),score);
    }
    
    public String toString() {
        return String.format("%s (%.3f)", getLocationString(), score);
    }
    
    public boolean equals(Object o) {
        if (o instanceof ScoredPoint) {
            ScoredPoint r = (ScoredPoint)o;
            if(!super.equals(r)) { return false; }
            if(score != r.score) { return false; }
            return true;
        } else {
            return false;
        }
    }
}
