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

public class DifferentialBindingGenerator implements Expander<Region,BindingExtent>, Closeable {

    private Expander<Region,BindingExtent> primaryGenerator;
    private Expander<Region,BindingExtent> backgroundGenerator;
    
    public DifferentialBindingGenerator(Expander<Region,BindingExtent> fg, 
            Expander<Region,BindingExtent> bg) { 
        primaryGenerator = fg;
        backgroundGenerator = bg;
    }
    
    public Iterator<BindingExtent> execute(Region a) {
        
        /*
         * Trying not to make the assumption that the extents come out of the underlying
         * generators sorted in any way.
         */
        
        Iterator<BindingExtent> primaryExtentItr = primaryGenerator.execute(a);
        Vector<BindingExtent> pextents = new Vector<BindingExtent>();
        while(primaryExtentItr.hasNext()) { pextents.add(primaryExtentItr.next()); }
                
        Iterator<BindingExtent> bgExtentItr = backgroundGenerator.execute(a);
        Vector<BindingExtent> bextents = new Vector<BindingExtent>();
        while(bgExtentItr.hasNext()) { pextents.add(bgExtentItr.next()); }
        
        BindingExtent[] bgarray = bextents.toArray(new BindingExtent[bextents.size()]);
        Arrays.sort(bgarray);
        
        Vector<BindingExtent> subtracted = new Vector<BindingExtent>();

        for(BindingExtent primary : pextents) { 
            subtracted.addAll(subtractExtents(primary, bgarray));
        }
        
        return subtracted.iterator();
    }
    
    private LinkedList<BindingExtent> subtractExtents(BindingExtent primary, BindingExtent[] bg) { 
        int i, j; 
        for(i = 0; i < bg.length && bg[i].getExtentEnd() < primary.getExtentStart(); i++) {}
        for(j = bg.length-1; j >= 0 && bg[j].getExtentStart() > primary.getExtentEnd(); j--) {}
        
        int minE = primary.getExtentStart(), maxE = primary.getExtentEnd();
        
        while(i < bg.length && bg[i].getExtentStart() <= primary.getExtentStart()) {
            minE = Math.max(minE, bg[i].getExtentEnd());
            i++;
        }
        
        while(j >= 0 && bg[j].getExtentEnd() >= primary.getExtentEnd()) { 
            maxE = Math.min(maxE, bg[j].getExtentStart());
            j--;
        }
        
        LinkedList<BindingExtent> fragments = new LinkedList<BindingExtent>();
        int k;
        for(k = i + 1; k < j; k++) { 
            BindingExtent ext = new BindingExtent(primary, minE, bg[k].getExtentStart()-1);
            minE = bg[k].getExtentEnd()+1;
            fragments.addLast(ext);
        }
        fragments.addLast(new BindingExtent(primary, minE, maxE));
        
        return fragments;
    }

    public void close() {
        if(primaryGenerator instanceof Closeable) { ((Closeable)primaryGenerator).close(); }
        if(backgroundGenerator instanceof Closeable) { ((Closeable)backgroundGenerator).close(); }
        primaryGenerator = backgroundGenerator = null;
    }

    public boolean isClosed() {
        return primaryGenerator == null && backgroundGenerator == null;
    }
    
}
