/**
 * 
 */
package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.warpdrive.model.*;

/**
 * @author tdanford
 */
public class SingleBindingEventPaintable extends RegionPaintable {
    private static ColorSet colors = new ColorSet();	
    private BindingEventModel model;
    private Hashtable<String,Object> typeColors;
    private int minWidth, minPixels;
    private SingleBindingEventProperties props;

    private Vector<BindingEvent> events;

    public SingleBindingEventPaintable(BindingEventModel m) {
        super();
        model = m;
        events = new Vector<BindingEvent>();
        minWidth = 1;
        minPixels = 8;
        props = new SingleBindingEventProperties();
        props.Color = colors.getColor(this.toString());
        model.addEventListener(this);
    }
    public SingleBindingEventProperties getProperties() {
        return props;
    }
    public void setLabel(String l) {
        props.Color = colors.getColor(l);
    }

    public BindingEventModel getModel() { return model; }
    public int getRecentTracks() { return 1; }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }
    
    public int getMaxVertSpace() { 
        int numTracks = 1;
        return Math.min(Math.max(20,numTracks * 8),60);        
    }
	
    public synchronized void eventRegistered(EventObject e) { 
        if (e.getSource() == model && model.isReady()) { 
            setCanPaint(true);
            setWantsPaint(true);
            
            Iterator<BindingEvent> evti = model.getResults();
            events.clear();
            while(evti.hasNext()) { events.add(evti.next()); }
            //layout.setRegions(events);
            
            notifyListeners();
        }		
    }
    
    public void paintItem(Graphics2D g, int x1, int y1, int x2, int y2) {
        if(!canPaint()) { return; }

        int w = x2 - x1, h = y2 - y1;
        clearLabels();

        //int[] xp = new int[3], yp = new int[3];

        int numTracks = 1;
        int track = 0;
        int trackHeight = numTracks > 0 ? Math.max(1, h / numTracks) : 1;
        int spacing = Math.max(2, trackHeight/10);
        Region region = model.getRegion();
        int start = region.getStart(), end = region.getEnd();

        // set the min-width.
        double basesPerPixel = (double)(end-start) / (double)w;
        minWidth = (int)Math.round(basesPerPixel * (double)minPixels);
        
        Vector<BindingWrapper> wrappers = new Vector<BindingWrapper>();
        for(BindingEvent evt : events) { wrappers.add(new BindingWrapper(evt, minWidth)); }

        for(int i = 0; i < wrappers.size(); i++) { 
            BindingWrapper wrapper = wrappers.get(i);

            BindingEvent evt = wrapper.event;
            String key = evt.getType();
            //if(!typeColors.containsKey(key)) { typeColors.put(key, colors.getColor()); }
            //color = (Color)typeColors.get(key);
			
            //int estart = evt.getStart(), eend = evt.getEnd();
            int estart = wrapper.getStart(), eend = wrapper.getEnd();
            
            int ex1 = getXPos(estart, start, end, x1, x2);
            int ex2 = getXPos(eend, start, end, x1, x2);
            
            int ewidth = Math.max(1, ex2-ex1);
            int ey1 = y1 + trackHeight * track;            
			
            g.setColor(props.Color);
            if(wrapper.expanded()) { 
                g.fillRect(ex1, ey1+spacing, ewidth, trackHeight-spacing*2);
                g.setColor(Color.black);
                g.drawRect(ex1, ey1+spacing, ewidth, trackHeight-spacing*2);
                addLabel(ex1, y1, ewidth, h, evt.getType());
            } else {
                g.fillRect(ex1, ey1+spacing, ewidth, trackHeight-spacing*2);
                g.setColor(Color.black);
                g.drawRect(ex1, ey1+spacing, ewidth, trackHeight-spacing*2);
                //g.drawRect(ex1, ey1, ewidth, trackHeight);
                addLabel(ex1, y1, ewidth, h, evt.getType());
                /*
                xp[0] = ex1; xp[1] = ex1 + (trackHeight/4); xp[2] = ex1 - (trackHeight/4);
                yp[0] = ey1; yp[1] = ey1 + trackHeight; yp[2] = ey1 + trackHeight;
                g.drawPolygon(xp, yp, 3);
                addLabel(ex1-(trackHeight/4), y1, trackHeight/2, h, evt.getType());
                */
            }
        }

    }
    
    public static class BindingWrapper extends Region { 
        
        public BindingEvent event;
        
        public BindingWrapper(BindingEvent evt, int minWidth) { 
            super(evt.getGenome(), evt.getChrom(), 
                    evt.getWidth() >= minWidth ? evt.getStart() : evt.getStart() - (minWidth-evt.getWidth())/2,
                    evt.getWidth() >= minWidth ? evt.getEnd() : evt.getEnd() + (minWidth-evt.getWidth())/2);
            event = evt;
        }
        
        public boolean expanded() { return getStart() < event.getStart() || getEnd() > event.getEnd(); }
        
        public int hashCode() { 
            int code = super.hashCode(); 
            code += event.hashCode(); code *= 37;
            return code;
        }
        
        public boolean equals(Object o) { 
            if(!(o instanceof BindingWrapper)) { return false; }
            BindingWrapper w = (BindingWrapper)o;
            if(!super.equals(w)) { return false; }
            if(!event.equals(w.event)) { return false; }
            return true;
        }
    }
}

