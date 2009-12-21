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
public class BayesDomainGenerator implements Expander<Region,BindingExtent>, Closeable {

    private BayesBindingGenerator bayes;
    private int minSize;
    private double density;

    public BayesDomainGenerator (BayesBindingGenerator b, int minSize, double density) {
        bayes = b;
        bayes.setPeaks(false);
        this.minSize = minSize;
        this.density = density;
    }

    public Iterator<BindingExtent> execute(Region region) {
        ArrayList<BindingExtent> output = new ArrayList<BindingExtent>();
        Iterator<BindingExtent> events = bayes.execute(region);
        ArrayList<Integer> positions = new ArrayList<Integer>();
        while (events.hasNext()) {
            Region r = events.next();
            positions.add((r.getStart() + r.getEnd()) / 2);
        }
        for (int i = 0; i < positions.size(); i++) {
            for (int j = positions.size() - 1; j > i; j--) {
            	int diff = positions.get(j) - positions.get(i);
                if (diff < minSize) {break;}
                
                //if ((j - i + 1) / (diff) < density) {continue;}
                if((double)diff / (double)(j-i+1) > density) { continue; }
                
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
                
                // this should be unnecessary, right?
                //break; 
            }
        }
        return output.iterator();
    }


    public boolean isClosed() {return bayes.isClosed();}
    public void close() {bayes.close();}
}

