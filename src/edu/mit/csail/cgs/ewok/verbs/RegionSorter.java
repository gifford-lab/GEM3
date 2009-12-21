/*
 * Created on Aug 9, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;

public class RegionSorter<X extends Region> implements Mapper<Iterator<X>,Iterator<X>> {
   
    public RegionSorter() {}

    public Iterator<X> execute(Iterator<X> a) {
        Vector<X> regions = new Vector<X>();
        while(a.hasNext()) { regions.add(a.next()); }
        Region[] array = regions.toArray(new Region[regions.size()]);
        Arrays.sort(array);
        return new RegionArrayIterator<X>(array);
    } 
}

class RegionArrayIterator<X extends Region> implements Iterator<X> {
    
    private Region[] array;
    private int index;
    
    public RegionArrayIterator(Region[] array) { 
        this.array = array;
        index = 0;
    }

    public boolean hasNext() {
        return array != null && index < array.length;
    }

    public X next() {
        X val = (X)array[index++];
        if(index >= array.length) { array = null; }
        return val;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    } 
}
