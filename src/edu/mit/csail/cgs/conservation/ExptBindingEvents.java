package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.Saveable;

/**
 * 
 * @author tdanford
 *
 * ExptBindingEvents is the class which contains the *raw* material for calling 
 * binding -- the total set of BindingEvent objects for a given ExptDescriptor.
 */

public class ExptBindingEvents implements Saveable {

	private ExptDescriptor expt;
	private Vector<BindingEvent> evts;
	
	public ExptBindingEvents(ExptDescriptor ed) { 
		expt = ed;
        evts = new Vector<BindingEvent>();
	}
	
	public ExptBindingEvents(Genome g, DataInputStream dis) throws IOException { 
		expt = new ExptDescriptor(dis);
        evts = new Vector<BindingEvent>();
        int s = dis.readInt();
        for(int i = 0; i < s; i++) { 
        	evts.add(new BindingEvent(g, dis));
        }
	}
	
	public void save(DataOutputStream dos) throws IOException {
		expt.save(dos);
		dos.writeInt(evts.size());
		for(BindingEvent evt : evts) { 
			evt.save(dos);
		}
	}
	
	public ExptBindingEvents(Genome g, ExptDescriptor e, PeakCaller pc) { 
        expt = e;
        evts = new Vector<BindingEvent>();

        Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
        Iterator<Region> chromRegions = 
        	new MapperIterator<NamedRegion,Region>(new CastingMapper<NamedRegion,Region>(), 
        			chroms);
        Iterator<BindingEvent> events = 
        	new ExpanderIterator<Region,BindingEvent>(pc, chromRegions);
        
        while(events.hasNext()) { addBindingEvent(events.next()); }
	}
    
    public ExptBindingEvents(Genome g, ExptDescriptor e, File f) throws IOException { 
        expt = e;
        evts = new Vector<BindingEvent>();
        
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line = null;
        while((line = br.readLine()) != null) { 
            line = line.trim();
            if(line.length() > 0 && line.charAt(0) != '#') { 
                String[] array = line.split("\\s+");
                String c = array[0];
                int strt = Integer.parseInt(array[1]), end = Integer.parseInt(array[2]);
                double str = Double.parseDouble(array[3]), conf = Double.parseDouble(array[4]);
                BindingEvent evt = new BindingEvent(g, c, strt, end, str, conf, e.getName());
                addBindingEvent(evt);
            }
        }
        br.close();
    }
	
	public void addBindingEvent(BindingEvent evt) { evts.add(evt); }
	public ExptDescriptor getExpt() { return expt; }
    
    public Vector<BindingEvent> getEventSubset(Region r) { 
        Vector<BindingEvent> subset = new Vector<BindingEvent>();
        for(BindingEvent e : evts) { if(e.overlaps(r)) { subset.add(e); } }
        return subset;
    }
	
	public Filter<Region,Region> getBindingFilter() { 
		return new BindingFilter();
	}
	
	private class BindingFilter implements Filter<Region,Region> {
		
		public BindingFilter() {}

		public Region execute(Region a) {
			for(BindingEvent evt : evts) { 
				if(a.overlaps(evt)) { 
					return a;
				}
			}
			return null;
		} 
		
	}
	
}
