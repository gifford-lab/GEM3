package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;

/**
 * A StaticExpander expands an X to all the objects in the input fixed
 * set of Ys that overlap the X.  The collection of Ys is stored in
 * memory and searched linearly when execute() is called.
 *
 * This class is useful when you have an in-memory collection of Regions (or subclasses)
 * that you wish to search for overlaps without getting them into the database.  It is not well suited to large
 * input sets.
 */
 
public class StaticExpander<X extends Region,Y extends Region> implements Expander<X,Y> {

    public ArrayList<Y> objects;

    public StaticExpander(Collection<Y> o) {
        objects = new ArrayList<Y>();
        objects.addAll(o);
        Collections.sort(objects);
    }
    public StaticExpander() {
        objects = new ArrayList<Y>();
    }
    public void add(Collection<Y> o) {
        objects.addAll(o);         
        Collections.sort(objects);
    }
    
    public Iterator<Y> execute(X region) {
        ArrayList<Y> out = new ArrayList<Y>();
        for (int i = 0; i < objects.size(); i++) {
            if (region.overlaps(objects.get(i))) {
                out.add(objects.get(i));
            }
        }
        return out.iterator();
    }

}