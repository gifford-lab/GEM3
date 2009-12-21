/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.Sequence;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 *
 * <code>ChromRegionIterator</code> is an iterator over the chromosomes
 * in the Genome provided to the constructor.
 */
public class ChromRegionIterator implements Iterator<NamedRegion> {
    
    private Genome genome;
    private Vector<NamedRegion> regions;
    private int index;

    public ChromRegionIterator(Genome g) {
        genome = g;
        regions = new Vector<NamedRegion>();
        index = 0;
        
        for(String chromName : genome.getChromList()) { 
            Genome.ChromosomeInfo s = genome.getChrom(chromName);
            NamedRegion region = new NamedRegion(genome, chromName, 0, s.getLength() - 1, chromName);
            regions.add(region);
        }
        
    }
    

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        return index < regions.size();
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public NamedRegion next() {
        NamedRegion nr = regions.get(index);
        index++;
        //System.out.println("--- " + nr.getChrom());
        return nr;
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

}
