package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.ScoredRegion;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;

// thin painter that draws ScoredRegions.  It draws them as rectangles and uses
// a shade of gray that's proportional to the score.

public class ShadedScoredPainter extends RegionPaintable {

    private double maxScore;
    private RegionExpanderModel<ScoredRegion> model;
    private PaintableProperties props;

    public ShadedScoredPainter(RegionExpanderModel<ScoredRegion> model, double maxScore) {
        super();
        this.model = model;
        model.addEventListener(this);
        this.maxScore = maxScore;
        props = new PaintableProperties();
    }
    public PaintableProperties getProperties() {
        return props;
    }
    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
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
        return 20;
    }
    
    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        if (!canPaint()) {
            return;
        }
        Iterator<ScoredRegion> regions = model.getResults();
        while (regions.hasNext()) {
            ScoredRegion r = regions.next();
            int minx = getXPos(r.getStart(),
                               getRegion().getStart(),
                               getRegion().getEnd(),
                               x1,
                               x2);
            int maxx = getXPos(r.getEnd(),
                               getRegion().getStart(),
                               getRegion().getEnd(),
                               x1,
                               x2);
            float percent = (float)(1.0 - r.getScore() / maxScore);
            g.setColor(new Color(percent,percent,percent,(float)1.0));
            g.fillRect(minx,y1,maxx-minx,(y2-y1));
        }
        if (props.DrawTrackLabel) {
            g.setColor(new Color((float)0.0,(float)0.0,(float)0.0,(float)0.5));
            g.drawString(getLabel(),x1,y2);
        }
    }
}
    
