/**
 * 
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import java.util.*;
import java.sql.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.binding.BindingScanLoader;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;


/**
 * @author Timothy Danford
 */
public class BindingExpander 
	implements Expander<Region,BindingEvent>, SelfDescribingVerb, Closeable {
	
	private boolean closeLoader;
    private Vector<BindingScan> scans;
    private BindingScanLoader loader;
    
    public BindingExpander() { 
    	try {
    		closeLoader = true;
			loader = new BindingScanLoader();
			scans = new Vector<BindingScan>();
			
		} catch (UnknownRoleException e) {
			e.printStackTrace();
			throw new IllegalStateException(e.getMessage());
		} catch (SQLException e) {
			e.printStackTrace();
			loader = null;
		}
    }
	
    public BindingExpander(BindingScanLoader l, BindingScan s) {
    	closeLoader = false;
        loader = l;
        scans = new Vector<BindingScan>();
        scans.add(s);
    }

    public BindingExpander(BindingScanLoader l, Collection<BindingScan> s) {
    	closeLoader = false;
        loader = l;
        scans = new Vector<BindingScan>(s);
    }
    
    public void close() { 
    	if(closeLoader) { 
    		loader.close();
    	}
    	loader = null;
    }
    
    public boolean isClosed() { 
    	return loader==null;
    }

    public void addBindingScan(BindingScan bs) { scans.add(bs); }

    public Iterator<BindingEvent> execute(Region a) {
        try { 
            LinkedList<BindingEvent> lst = new LinkedList<BindingEvent>();
			
            for(BindingScan scan : scans) { 
                Collection<BindingEvent> events = loader.loadEvents(scan, a);
                lst.addAll(events);
            }
            
            return lst.iterator();
        } catch(SQLException se) { 
            se.printStackTrace(System.err);
            return new EmptyIterator<BindingEvent>();
        }
    }
    
    private static final String[] inputNames = { "Regions" };
    private static final EchoType[] inputClasses = { new ClassType(Region.class) };
    private static final String[] paramNames = { "Scan" };
    private static final EchoType[] paramClasses = { new ClassType(BindingScan.class) };
    
	public EchoType[] getInputClasses() { return inputClasses; }
	public String[] getInputNames() {return inputNames; }

	public EchoType getOutputClass() {
		return new ClassType(BindingEvent.class);
	}

	public EchoType[] getParameterClasses() { return paramClasses; }
	public String[] getParameterNames() { return paramNames; }

	public void init(Map<String, Object> params) {
		BindingScan bs = (BindingScan)params.get(paramNames[0]);
		scans.clear();
		scans.add(bs);
	}
}
