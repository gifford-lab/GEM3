package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.model.RegionMapperModel;


public class PerBaseScorePainter<X extends Number> extends RegionPaintable {


    private RegionMapperModel<X[]> model;
    private PerBaseScoreProperties props;

    /* the assumption is that maxscore > 0 and minscore < 0 */
    public PerBaseScorePainter(RegionMapperModel<X[]> model,
                               double minscore,
                               double maxscore) {
        this.model = model;
        props = new PerBaseScoreProperties();
        props.MaxScore = maxscore;
        props.MinScore = minscore;
        props.MiddleScore = 0.9;
        model.addEventListener(this);
    }

    public PerBaseScorePainter(RegionMapperModel<X[]> model,
                               double minscore,
                               double middlescore,
                               double maxscore) {
        this.model = model;
        props = new PerBaseScoreProperties();
        props.MaxScore = maxscore;
        props.MinScore = minscore;
        props.MiddleScore = middlescore;
        model.addEventListener(this);
    }
    public PerBaseScoreProperties getProperties () {return props;}

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
        g.setColor(Color.WHITE);
        g.fillRect(x1,y1,x2-x1,y2-y1);
        if (!canPaint()) {
            return;
        }
        boolean usemax = props.UseMax;
        int w = x2 - x1;
        double h = y2 - y1;
        Region region = getRegion();
        int regionstart = region.getStart(), regionend = region.getEnd();
        int regionwidth = regionend - regionstart;
        int pixwidth = props.PixWidth;
        if (pixwidth < 1) {
            pixwidth = 1;
        }
        X[] vals = model.getResults();
        for (int i = x1; i < x2; i += pixwidth) {
            int rstart = (int)Math.round((i - x1) * regionwidth / ((double)w));
            int rend = (int)Math.round((i - x1 + pixwidth) * regionwidth / ((double)w));
            if (rstart < 0) { rstart = 0;}
            if (rend >= vals.length) {rend = vals.length - 1;}

            double sum = 0;
            if (usemax) {
                sum = Double.NEGATIVE_INFINITY;
                for (int j = rstart; j <= rend ; j++) {
                    sum = Math.max(vals[j].doubleValue(),sum);
                }
            } else {
                for (int j = rstart; j <= rend ; j++) {
                    sum += vals[j].doubleValue();
                }
                sum /= rend - rstart + 1;
            }            
            if (sum > props.MiddleScore) {
                g.setColor(Color.BLACK);
                sum = (sum - props.MiddleScore) / (props.MaxScore - props.MiddleScore);
            } else {
                g.setColor(Color.RED);
                sum = (props.MiddleScore - sum) / (props.MiddleScore - props.MinScore);
            }
            int fill = (int) (sum * h);
            g.fillRect(i,y1,pixwidth,fill);
        }
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLUE);
            g.drawString(getLabel(),x1,y2);
        }
    }
}

