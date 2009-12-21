package edu.mit.csail.cgs.warpdrive.paintable;

import java.util.EventObject;
import java.util.Set;
import java.util.ArrayList;
import java.awt.*;

import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.warpdrive.model.SpeciesAlignModel;
import edu.mit.csail.cgs.utils.*;

/* this painter goes with SpeciesAlignPainter (and the code is very similar) 
   but draws inside the normal hyperdrive window instead of the separate
   alignment window.  It draws colored bars to show where the alignment regions
   from the alignment window fall in the main window */
public class SpeciesAlignThinPainter extends RegionPaintable implements Listener<EventObject> {    

    private SpeciesAlignModel model;
    private DynamicAttribute attrib;
    private Font myFont;
    private ColorSet cs;
    private HashMarkPaintable hmp;
    private int lastpoint = 0;
    private Genome myGenome;
    private PaintableProperties props;

    public SpeciesAlignThinPainter(SpeciesAlignModel model, Genome g) {
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        myFont = new Font("Arial", Font.PLAIN, 12);
        cs = new ColorSet();
        hmp = new HashMarkPaintable();
        props = new PaintableProperties();
        myGenome = g;
    }
    public PaintableProperties getProperties() {return props;}

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public int getMaxVertSpace() {return 30;}

    public void eventRegistered(EventObject e) {
        if (e.getSource() == model &&
            model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        if (!model.isReady()) {
            return;
        }
        Set<Genome> genomes = model.getGenomes();
        if (genomes == null) {
            g.drawString("No Genomes Yet",x1+40,y1+40);
            return;
        }
        Genome currentGenome = model.getCurrentGenome();
        if (currentGenome == null) {
            g.drawString("Haven't moved",x1+40,y1+40);
            return;
        }
        genomes.remove(currentGenome);
        if (genomes.size() == 0) {
            g.drawString("Single genome isn't very interesting",x1+40,y1+40);
            return;
        }
        cs.reset();
        int h = y2 - y1;
        int midY = y1 + h/2;
        int rstart = getRegion().getStart();
        int rend = getRegion().getEnd();
        Stroke old = g.getStroke();
        Stroke thin = new BasicStroke(1);
        Stroke medium = new BasicStroke(3);
        Stroke thick = new BasicStroke(6);
        for (Genome genome : genomes) {
            Region best = model.getBestRegion(genome);
            if (best == null) {continue;}
            String bestc = best.getChrom();
            int ostart = best.getStart();
            int oend = best.getEnd();
            java.util.List<MultiZAlignRegion> alignedRegions = model.getAlignedRegions(genome);
            for (MultiZAlignRegion align : alignedRegions) {
                if (align.getOtherChrom().equals(bestc)) {
                    g.setColor(cs.colorAt(hashPosPair(align.getStart(),
                                                      align.getEnd(),
                                                      align.getOtherStart(),
                                                      align.getOtherEnd(),
                                                      cs.colorCount())));
                    int rsx = getXPos(align.getStart(),rstart,rend,x1,x2);
                    int osx = getXPos(align.getOtherStart(),ostart,oend,x1,x2);
                    int rex = getXPos(align.getEnd(),rstart,rend,x1,x2);
                    int oex = getXPos(align.getOtherEnd(),ostart,oend,x1,x2);


                    g.setStroke(medium);
                    if (myGenome == currentGenome) {
                        g.drawLine(rsx,midY,rex,midY);                        
                        g.drawLine(rsx,y1,rsx,y2);
                        g.drawLine(rex,y1,rex,y2);
                        addLabel(rsx, y1, rex - rsx, y2, align.toString());
                    } else {
                        g.drawLine(osx,midY,oex,midY);                        
                        g.drawLine(osx,y1,osx,y2);
                        g.drawLine(oex,y1,oex,y2);
                        addLabel(osx, y1, oex - osx, y2, align.toString());
                    }

                }
            }
        }
        g.setStroke(old);
        g.setColor(Color.BLACK);
        g.drawString(model.getAlignment(),x1+5,y2);
    }

    public int hashPosPair (int x1, int x2, int x3, int x4, int max) {
        return Math.abs(x1 + x2 * 5 + x3 * 3 + x4 * 11) % max;
    }

}
