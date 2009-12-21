/*
 * Created on Apr 4, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Vector;

import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.Closeable;

public class CollapsingBindingGenerator implements Expander<Region,BindingExtent>, Closeable {
    
    private Expander<Region,BindingExtent> gen;
    
    public CollapsingBindingGenerator(Expander<Region,BindingExtent> g) { 
        gen = g;
    }

    public Iterator<BindingExtent> execute(Region a) {
        Vector<BindingExtent> extents = new Vector<BindingExtent>();
        Iterator<BindingExtent> itr = gen.execute(a);
        while(itr.hasNext()) { extents.add(itr.next()); }
        
        LinkedList<BindingExtent> newextents = new LinkedList<BindingExtent>();
        
        BindingExtent[] array = extents.toArray(new BindingExtent[extents.size()]);
        Arrays.sort(array);
        
        for(int i = 0; i < array.length; i++) {
            int start = i;
            int end = i+1;
            while(end < array.length && array[end].overlaps(array[start])) { 
                end++;
            }
            
            int estart = array[start].getExtentStart();
            int eend = array[end-1].getExtentEnd();
            
            BindingExtent ext = new BindingExtent(array[start], estart, eend);
            newextents.addLast(ext);
            
            i = end-1;
        }
        
        return newextents.iterator();
    }

    public void close() {
        if(gen instanceof Closeable) { 
            ((Closeable)gen).close();
        }
        gen = null;
    }

    public boolean isClosed() {
        return gen == null;
    }
    
}
