package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.Closeable;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Iterator;

/* Domains of binding are those regions in which the density of underlying binding events (bayes binding events, in this case)
   exceed some threshold and the size of the domain exceeds a minimum size.
*/
public class BindingDomainGenerator implements Expander<Region,BindingExtent>, Closeable {

    private Expander<Region,BindingExtent> generator;
    private int minSize;
    private double density;

    public BindingDomainGenerator(Expander<Region,BindingExtent> b, int minSize, double density) {
        generator = b;
        this.minSize = minSize;
        this.density = density;
    }

    public Iterator<BindingExtent> execute(Region region) {
        ArrayList<BindingExtent> output = new ArrayList<BindingExtent>();
        Iterator<BindingExtent> events = generator.execute(region);
        ArrayList<Integer> positions = new ArrayList<Integer>();
        while (events.hasNext()) {
            Region r = events.next();
            positions.add((r.getStart() + r.getEnd()) / 2);
        }
        for (int i = 0; i < positions.size(); i++) {
            for (int j = positions.size() - 1; j > i; j--) {

                int posDiff = positions.get(j)-positions.get(i);
                int indexDiff = j-i+1;
                
                if (posDiff < minSize) {break;}
                
                //if ((j - i + 1) / (positions.get(j) - positions.get(i)) < density) {continue;}
                if((double)posDiff / (double)indexDiff > density) { continue; }
                
                output.add(new BindingExtent(region.getGenome(),
                                             region.getChrom(),
                                             positions.get(i),
                                             positions.get(j),
                                             positions.get(j) - positions.get(i),
                                             (j - i + 1) / (positions.get(j) - positions.get(i)),
                                             "Domain",
                                             positions.get(i),
                                             positions.get(j)));
                i = j + 1;
                
                // shouldn't be needed, ince j = (i-1), but then +1 at the 
                // end of the loop means j=i, which will fail the continuation test.
                //break;
            }
        }
        return output.iterator();
    }

    public void close() {
        if(generator instanceof Closeable) { 
            ((Closeable)generator).close();
            generator = null;
        }
    }

    public boolean isClosed() {
        return generator==null;
    }
}

