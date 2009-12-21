package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionMapperModel;


/* paints some subset of nucleotides blue and some subset red.  Default
   is blue = GC and red = AT
*/
public class GCContentPainter extends RegionPaintable {

    private RegionMapperModel<String> model;
    private GCContentProperties props;
    
    public GCContentPainter(RegionMapperModel<String> model) {
        super();
        this.model = model;
        props = new GCContentProperties();
        model.addEventListener(this);
    }
    public GCContentProperties getProperties() { return props;}
    
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
            //            System.err.println("GCCP got event and notified");
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
        int w = x2 - x1;
        double h = y2 - y1;
        Region region = getRegion();
        int regionstart = region.getStart(), regionend = region.getEnd();
        int regionwidth = regionend - regionstart;
        int pixwidth = props.PixWidth;
        if (pixwidth < 1) {
            pixwidth = 1;
        }

        String wholestring = model.getResults();
        if (wholestring == null) {
            return;
        }
        char[] chars = wholestring.toUpperCase().toCharArray();
        String bluestring = props.BlueBases;
        char[] bluechar = new char[bluestring.length()];        
        for (int i = 0; i < bluechar.length; i++) {
            bluechar[i] = bluestring.charAt(i);
        }
        String redstring = props.RedBases;
        char[] redchar = new char[redstring.length()];
        
        for (int i = 0; i < redchar.length; i++) {
            redchar[i] = redstring.charAt(i);
        }

        float expfrac = (float)bluechar.length  / ((float)redchar.length + (float)bluechar.length);
        
        for (int i = x1; i < x2; i += pixwidth) {
            int rstart = (int)Math.round((i - x1) * regionwidth / ((double)w));
            int rend = (int)Math.round((i - x1 + pixwidth) * regionwidth / ((double)w));
            if (rend >= chars.length) {
                rend = chars.length - 1 ;
            }
            int blue = 0; int red = 0;
            for (int j = rstart; j <= rend; j++) {
                for (int k = 0; k < bluechar.length; k++) {
                    if (chars[j] == bluechar[k]) {
                        blue++;
                        break;
                    }
                }
                for (int k = 0; k < redchar.length; k++) {
                    if (chars[j] == redchar[k]) {
                        red++;
                        break;
                    }
                }
            }
            if (red + blue == 0) {continue;}               
            float frac = (float)blue / ((float)red + (float) blue);
            if (frac > expfrac) {
                g.setColor(Color.BLUE);
                int fill = (int) (h * ((frac - expfrac) / (1.0 - expfrac)));
                g.fillRect(i,y1,pixwidth,fill);
            } else {
                g.setColor(Color.RED);
                int fill = (int)(h * ((expfrac - frac) / expfrac));
                g.fillRect(i,y1,pixwidth,fill);
            }
        }
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1,y2);
        }
    }
}

