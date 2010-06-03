package edu.mit.csail.cgs.deepseq.features;

import java.util.ArrayList;
import java.util.Collection;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.Gene;

/**
 * Any genomic feature that we wish to discover in a deep-sequencing experiment. 
 * 
 * @author shaun
 *
 */
public abstract class Feature{
	public Region coords;
	public Gene nearestGene=null;
	public int distToGene=Integer.MAX_VALUE;
    public Collection<Region> annotations;
	
    public Feature(){this(null);}
	public Feature(Region c){
		coords=c;
	}
	
	public Feature(Region c, Gene g, int dist, Collection<Region> annots) {
		coords = c; nearestGene = g; distToGene = dist; annotations = annots;
	}
	
	public Feature(Region c, Gene g, int dist) {
		coords = c; nearestGene = g; distToGene = dist;
	}
	
    public void addAnnotation(Region r) {
        if (annotations == null) {
            annotations = new ArrayList<Region>();
        }
        annotations.add(r);
    }	
    
    public Point getPeak(){return new Point(coords.getGenome(), coords.getChrom(), (coords.getStart()+coords.getEnd())/2);}

    public abstract String toString();
    public abstract String toGFF();
    public abstract String headString();
    public abstract String toSequence(int win);
}
