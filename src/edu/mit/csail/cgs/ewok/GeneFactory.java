/*
 * Created on Sep 28, 2006
 */
package edu.mit.csail.cgs.ewok;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.Expander;

/**
 * @author tdanford
 *
 * <code>GeneFactory</code> returns an Expander than maps Regions to Genes.  The purpose
 * of the Factory is to return an appropriate expander for a given Genome.
 */
public interface GeneFactory {
    public Expander<Region,Gene> getExpander(Genome g);
}
