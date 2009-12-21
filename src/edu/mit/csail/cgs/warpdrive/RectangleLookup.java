/*
 * Created on Nov 15, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.Pair;

public class RectangleLookup<X> {

    private Vector<Pair<X,Rectangle>> values;
    private Map<Rectangle,Set<X>> lookup;
    
    public RectangleLookup() { 
        values = new Vector<Pair<X,Rectangle>>();
        lookup = new HashMap<Rectangle,Set<X>>();
    }
    
    public void clear() { values.clear(); lookup.clear(); }
    
    public Collection<Pair<X,Rectangle>> getValueList() { 
        return new LinkedList<Pair<X,Rectangle>>(values); 
    }
    
    public void clearRectangle(Rectangle r) { 
        if(lookup.containsKey(r)) { 
            lookup.remove(r);
        }

        Iterator<Pair<X,Rectangle>> itr = values.iterator();
        while(itr.hasNext()) { 
            Pair<X,Rectangle> p = itr.next();
            if(p.getLast().equals(r)) { itr.remove(); }
        }
    }
    
    public void addValue(X v, Rectangle r) { 
        values.add(new Pair<X,Rectangle>(v, r));
        if(!lookup.containsKey(r)) { lookup.put(r, new HashSet<X>()); }
        lookup.get(r).add(v);
    }
    
    public Collection<X> getRectangleValues(Rectangle r) { 
        return lookup.get(r);
    }
    
    public X getFirstValue(Point p) { 
        for(Pair<X,Rectangle> v : values) { 
            if(v.getLast().contains(p)) { 
                return v.getFirst();
            }
        }
        return null;
    }
    
    public Collection<X> getAllValues(Point p) { 
        LinkedList<X> ret = new LinkedList<X>();
        for(Pair<X,Rectangle> v : values) { 
            if(v.getLast().contains(p)) { 
                ret.addLast(v.getFirst());
            }
        }
        return ret;
    }
}
