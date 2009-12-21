/*
 * Created on Dec 4, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.*;

import javax.swing.JProgressBar;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.*;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;

/* annotates a genome-wide set of binding events (provided by the Expander<Region,BindingEvent>)
   with refgenes
*/

public class BindingEventAnnotations extends CachedAnnotations<BindingEvent,Region> {
    
    private Expander<BindingEvent,? extends Region> annotator;
    private JProgressBar progressBar;
    
    public BindingEventAnnotations(Genome g, Expander<Region,BindingEvent> evts) {
        super();
        progressBar = null;
        annotator = null;
        init(evts,g);
    }

    /* annotator takes a binding event and returns the relevant annotations (genes, etc)
       evts returns the BindingEvents in a region; it is applied
         to all of the chromosomes to get the set of BindingEvents to be annotated
    */
    public BindingEventAnnotations(Genome g, 
    		Expander<BindingEvent,Region> annotator, 
    		Expander<Region,BindingEvent> evts, 
    		JProgressBar pbar) { 
        super();
        progressBar = pbar;
        this.annotator = annotator;
        init(evts,g);
    }
    
    public void init(Expander<Region,BindingEvent> evts, Genome g) {

        ChromRegionIterator chroms = new ChromRegionIterator(g);
        Iterator<Region> regionChroms = 
            new MapperIterator<NamedRegion,Region>(new CastingMapper<NamedRegion,Region>(), chroms);
        Iterator<BindingEvent> evtItr = new ExpanderIterator<Region,BindingEvent>(evts, regionChroms);
        
        LinkedList<BindingEvent> totalEvents = new LinkedList<BindingEvent>();
        while(evtItr.hasNext()) { 
            totalEvents.addLast(evtItr.next());
        }
        
        addItems(totalEvents);
        
        if(progressBar != null) { 
            progressBar.setMaximum(getNumItems());
        }

        if (annotator == null) {
            annotator = new RefGeneGenerator<BindingEvent>(g);
        }
        
        addAnnotations("genes", annotator);
    }
    
    public void markProgress() { 
        if(progressBar != null) {
            progressBar.setValue(progressBar.getValue() + 1);
        }
    }
    
    public Expander<? extends Region,? extends Region> getGeneGenerator() { return annotator; }
}
