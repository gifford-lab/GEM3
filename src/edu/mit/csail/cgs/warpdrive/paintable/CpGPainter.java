package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionMapperModel;

public class CpGPainter extends RegionPaintable {

    private RegionMapperModel<String> model;
    private CpGProperties props;
    
    public CpGPainter(RegionMapperModel<String> model) {
        super();
        this.model = model;
        props = new CpGProperties();
        model.addEventListener(this);
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public CpGProperties getProperties() {
        return props;
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
        //        System.err.println("chars.length = " + chars.length + "  regionwidth=" + regionwidth);
        g.setColor(Color.BLUE);
        for (int i = x1; i < x2; i += pixwidth) {
            int rstart = (int)Math.round((i - x1) * regionwidth / ((double)w));
            int rend = (int)Math.round((i - x1 + pixwidth) * regionwidth / ((double)w));
            if (rend >= chars.length) {
                rend = chars.length - 1 ;
            }
            int found = 0;
            for (int j = rstart; j < rend; j++) {
                if ((chars[j] == 'C' && chars[j+1] == 'G')) {
                    found++;
                }
            }
            float frac = 2 * ((float)found) / (rend - rstart);
            if (frac > 1) {
                frac = 1;
            }

            int fill = (int) (h * frac);
            //            System.err.println(rstart  + " - " + rend + "  : " + found + ", " + frac + ", " + fill);
            g.fillRect(i,y1,pixwidth,fill);
        }
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1,y2);
        }
    }
}

