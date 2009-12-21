/*
 * Created on May 21, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs.chippet;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chippet.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class ChipPetExpander implements Expander<Region,ChipPetDatum>, Closeable {
    
    private ChipPetLoader loader;
    private ChipPetExpt expt;
    private boolean shouldCloseLoader;
    
    public ChipPetExpander(ChipPetLoader loader, ChipPetExpt expt) { 
        this.loader = loader;
        this.expt = expt;
        shouldCloseLoader = false;
    }
    
    public ChipPetExpander(String exptName) throws SQLException {
        loader = new ChipPetLoader();
        expt = loader.loadExperiment(exptName);
        shouldCloseLoader = true;
    }

    public Iterator<ChipPetDatum> execute(Region a) {
        try {
            Collection<ChipPetDatum> intvs = loader.loadIntervals(expt, a);
            return intvs.iterator();
        } catch (SQLException e) {
            e.printStackTrace();
            return new EmptyIterator<ChipPetDatum>();
        }
    }

    public void close() {
        if(shouldCloseLoader) { 
            loader.close();
        }
        expt = null;
        loader = null;
    }

    public boolean isClosed() {
        return loader == null || loader.isClosed();
    }
    
    public Expander<Region,BindingEvent> getBindingEventWrapper(int thresh) { 
        return new BindingEventWrapper(thresh);
    }

    private class BindingEventWrapper implements Expander<Region,BindingEvent>, Closeable {
        
        private int threshold;
        
        public BindingEventWrapper(int t) { 
            threshold = t;
        }

        public void close() {
            ChipPetExpander.this.close();
        }

        public boolean isClosed() {
            return ChipPetExpander.this.isClosed();
        }

        public Iterator<BindingEvent> execute(Region a) {
            Iterator<ChipPetDatum> itr = ChipPetExpander.this.execute(a);
            LinkedList<BindingEvent> evts = new LinkedList<BindingEvent>();
            
            String type = "ChipPetInterval";
            double conf = 1.0;
            double size = (double)threshold;
            
            RunningOverlapSum rs = new RunningOverlapSum(a.getGenome(), a.getChrom());
            
            while(itr.hasNext()) { 
                ChipPetDatum intv = itr.next();
                for(int i = 0; i < intv.getCount(); i++) { 
                    rs.addRegion(intv);
                }
            }
                
            for(Region r : rs.collectRegions(threshold)) { 
                BindingEvent evt = new BindingEvent(a.getGenome(), a.getChrom(), 
                        r.getStart(), r.getEnd(), size, conf, type);
                evts.addLast(evt);
            }
            
            return evts.iterator();
        } 
        
    }
}
