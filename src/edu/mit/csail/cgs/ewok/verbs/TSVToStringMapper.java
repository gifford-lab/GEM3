package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;

/**
 * <code>TSVToStringMapper</code> tries to hide (or at least put all
 * the messy parts in one place) the complexity of mapping objects to
 * tab-separated-values as one might want for file output.  This works
 * by having a hardcoded list of classes and a block of code to output
 * each.
 */
public class TSVToStringMapper implements Mapper<Object, String> {

    public String execute(Object o) {
        if (o instanceof Genome) {
            Genome g = (Genome) o;
            return g.getSpecies() + "\t" + g.getVersion();
        } else if (o instanceof Gene) {
            Gene g = (Gene) o;
            StringBuffer out = new StringBuffer();
            out.append(g.getID() + "\t" + g.getName());
            for (String a : g.getAliases()) {
                out.append("\t" + a);
            }
            return out.toString();
        } else if (o instanceof Region) {
            Region r = (Region) o;
            return r.getChrom() + "\t" + r.getStart() + "\t" + r.getEnd();
        } else {
            return o.toString();
        }

    }

}