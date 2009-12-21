/*
 * Created on Dec 5, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.Region;

import java.sql.SQLException;
import java.util.*;

public class BindingScanGenerator implements Expander<Region,BindingEvent> {

    private BindingScanLoader loader;
    private LinkedList<BindingScan> scans;
    
    public BindingScanGenerator(BindingScanLoader l, BindingScan s) { 
        loader = l;
        scans = new LinkedList<BindingScan>(); scans.add(s);
    }

    public BindingScanGenerator(BindingScanLoader l) { 
        loader = l;
        scans = new LinkedList<BindingScan>(); 
    }
    
    public void addBindingScan(BindingScan bs) { scans.addLast(bs); }

    public Iterator<BindingEvent> execute(Region a) {
        LinkedList<BindingEvent> evts = new LinkedList<BindingEvent>();
        
        for(BindingScan scan : scans) { 
            try {
                Collection<BindingEvent> events = loader.loadEvents(scan, a);
                evts.addAll(events);
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
        
        return evts.iterator();
    }
}
