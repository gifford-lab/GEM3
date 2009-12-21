package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.ScoredStrandedRegion;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;

public class MotifScanPainter extends RegionPaintable {
    
    private RegionExpanderModel<ScoredStrandedRegion> model;
    private MotifScanProperties props;

    public MotifScanPainter (RegionExpanderModel<ScoredStrandedRegion> model) {
        super();
        this.model = model;
        props = new MotifScanProperties();
        model.addEventListener(this);
    }
    public MotifScanProperties getProperties() {return props;}

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

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }


    public int getMaxVertSpace() {return 25;}

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        if (!canPaint()) {
            return;
        } 
        int height = Math.abs(lry - uly);
        int half = height / 2;
        int center = uly + half;
        int minPixelWidth = props.MinPixelWidth;
        g.setColor(Color.BLACK);
        g.drawLine(ulx,center,lrx,center);
        g.setColor(new Color(107,66,38)); // brown.  don't use red because that's not a match in the per-base scanner
        double heightmult = props.HeightMultiplier;
        Iterator<ScoredStrandedRegion> iter = model.getResults();
        while (iter.hasNext()) {
            ScoredStrandedRegion r = iter.next();
            int x1 = getXPos(r.getStart(),
                             model.getRegion().getStart(),
                             model.getRegion().getEnd(),
                             ulx,lrx);
            int x2 = getXPos(r.getEnd(),
                             model.getRegion().getStart(),
                             model.getRegion().getEnd(),
                             ulx,lrx);
            if (x2 - x1 < minPixelWidth) {
                int missing = minPixelWidth - (x2 - x1);
                int h = (missing +1) /2;
                x2 += h;
                x1 -= h;
            }
            int h = (int)(r.getScore() * heightmult);
            if (h > half) {
                h = half;
            }
            if (r.getStrand() == '+') {
                g.fillRect(x1,center-h,x2-x1,h);
            } else {
                g.fillRect(x1,center,x2-x1,h);
            }
        }
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            Font oldFont = g.getFont();
            DynamicAttribute attrib = DynamicAttribute.getGlobalAttributes();
            Font nf = attrib.getLargeLabelFont(lrx - ulx,lry - uly);
            g.setFont(nf);
            g.drawString(getLabel(),ulx,lry);
            g.setFont(oldFont);
        }
    }

    
}

