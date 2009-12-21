package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;

public class PeakPainter extends RegionPaintable {
    
    private RegionExpanderModel<BindingEvent> model;
    private NonOverlappingLayout<BindingEvent> layout;
    private Font myFont;
    private DynamicAttribute attrib;
    private double htRat, wdRat;
    private PaintableProperties props;
    
    public PeakPainter(RegionExpanderModel<BindingEvent> model) {
        super();
        layout = new NonOverlappingLayout<BindingEvent>();
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        myFont = new Font("Arial", Font.PLAIN, 12);
        htRat = .03;
        wdRat = .03;
        props = new PaintableProperties();
        model.addEventListener(this);
    }
    public PaintableProperties getProperties() {
        return props;
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public synchronized void eventRegistered(EventObject e) {
        if (e.getSource() == model &&
            model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }

    public int getMaxVertSpace() { 
        return 40; 
    }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        if (!canPaint()) {
            return;
        }
        
        Iterator<BindingEvent> itr = model.getResults();
        Vector<BindingEvent> peaks = new Vector<BindingEvent>();
        while(itr.hasNext()) { peaks.add(itr.next()); }
        layout.setRegions(peaks);
        
        int w = x2 - x1, h = y2 - y1;
        int my = y1 + (h / 2);
 
        int numTracks = layout.getNumTracks();
        //System.out.println("Num Tracks: " + numTracks);
        
        int trackHeight = h;

        if(numTracks > 1) { 
            trackHeight = (int)Math.floor((double)h / (double)numTracks);
        }
        
        int geneHeight = Math.max(2, (int)Math.floor((double)trackHeight * 0.80));
        int halfHeight = trackHeight / 2;
        int halfGeneHeight = geneHeight / 2;

        for(int track=0; track<numTracks; track++) {
            g.setColor(Color.black);
            my = y1 + (halfHeight + (track * trackHeight));
            g.drawLine(x1, my, x2, my);
        }
        
        Region region = model.getRegion();
        int rs = region.getStart(), re = region.getEnd();
        int rw = re - rs + 1;

        // Additions by Lindy
        // Dynamically change font size as the size of the window changes
        int fontSize = attrib.getFontSize(w,h);
        fontSize = Math.max(7, (int)Math.floor((double)fontSize * 0.70));
        Font oldFont = g.getFont();
        g.setFont(new Font(myFont.getFontName(), myFont.getStyle(), fontSize));
        //--------------------------------------------------
        
        double xScale = (double)w / (double)rw;
        
        for(BindingEvent peak : peaks) { 
            int gs = peak.getStart(), ge = peak.getEnd();
            int track = 0;
            if(!layout.hasTrack(peak)) { 
                System.err.println("No track assigned to peak: " + peak.toString());
            } else { 
                track = layout.getTrack(peak);
            }
            
            my = y1 + (halfHeight + (track * trackHeight));
            
            int gx1 = x1 + (int)Math.round((double)(gs - rs) * xScale);
            int gx2 = x1 + (int)Math.round((double)(ge - rs) * xScale);
            int gxw = gx2 - gx1;
            int depth = Math.min(halfGeneHeight, gxw/3);

            int gleft = Math.max(gx1, x1), gright = Math.min(gx2, x2);
            int gtop = my - (halfGeneHeight/2), gbottom = gtop + (geneHeight/2);

            //g.setColor(Color.white);
            //g.fillRect(gx1, my - (halfGeneHeight / 2), gxw, geneHeight / 2);
            g.setColor(Color.pink);
            int gh = gbottom - gtop, gw = gright - gleft;
            int arc = Math.max(1, gh/2);
            g.fillRoundRect(gleft, gtop, gw, gh, arc, arc);  
            
            g.setColor(Color.black);
            //g.drawRect(gx1, my - (halfGeneHeight / 2), gxw, geneHeight / 2);
            g.drawLine(gleft, gtop, gright, gtop);
            g.drawLine(gleft, gbottom, gright, gbottom);
            if(gx1 >= x1) { g.drawLine(gx1, gtop, gx1, gbottom); }
            if(gx2 <= x2) { g.drawLine(gx2, gtop, gx2, gbottom); }
            
            // Somewhat prettier gene-name output
            //g.drawString(gene.getName(), gx1 + (gxw/3), my);
            int nx = Math.max(x1 + 3, gx1 + 3);  // gotta do the Math.max(), to make sure the name is on the screen.
            int ny = my + (halfGeneHeight / 2) - 2;
            //g.drawString(gene.getName(), nx, ny);
        }

        g.setFont(oldFont);
    }

    
}
