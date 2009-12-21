package edu.mit.csail.cgs.warpdrive.paintable;

import java.util.*;
import java.awt.*;

import edu.mit.csail.cgs.datasets.alignments.MultiZAlignRegion;
import edu.mit.csail.cgs.datasets.alignments.AlignmentStitcher;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.warpdrive.model.SpeciesAlignModel;
import edu.mit.csail.cgs.utils.*;


public class SpeciesAlignPainter extends RegionPaintable implements Listener<EventObject> {    

    private SpeciesAlignModel model;
    private DynamicAttribute attrib;
    private Font myFont;
    private ColorSet cs;
    private HashMarkPaintable hmp;
    private int lastpoint = 0;
    private SpeciesAlignProperties props;

    public SpeciesAlignPainter(SpeciesAlignModel model) {
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        myFont = new Font("Arial", Font.PLAIN, 12);
        cs = new ColorSet();
        hmp = new HashMarkPaintable();   
        props = new SpeciesAlignProperties();
    }
    public SpeciesAlignProperties getProperties() {return props;}

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

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
        if (props.PaintAligned) {
            paintAligned(genomes,currentGenome,g,x1,y1,x2,y2);
        } else {
            paintMismatch(genomes,currentGenome,g,x1,y1,x2,y2);
        }
    }

    private void paintAligned(Set<Genome> genomes,
                              Genome currentGenome,
                              Graphics2D g,
                              int x1, int y1, 
                              int x2, int y2) {
        cs.reset();
        int pixelsPerTrack = (y2 -  y1) / (1 + genomes.size());
        int rstart = getRegion().getStart();
        int rend = getRegion().getEnd();
        Font oldFont = g.getFont();
        int fontSize = attrib.getFontSize(x2-x1,y2-y1);
        fontSize = Math.max(8, (int)Math.floor((double)fontSize));
        g.setFont(new Font(myFont.getFontName(), myFont.getStyle(), fontSize));
        g.setColor(Color.BLACK);
        int baseY = y2 - pixelsPerTrack/2;
        // this is the base 
        Stroke oldStroke = g.getStroke();
        Stroke thin = new BasicStroke(1);
        Stroke medium = new BasicStroke(3);
        Stroke thick = new BasicStroke(6);
        int hmWidth = 20;
        int horizLineOffset = 6;

        hmp.setLabelAbove(false);
        hmp.setRegion(getRegion());
        hmp.paintItem(g,x1,baseY,x2,baseY+hmWidth);
        g.drawString(currentGenome.getVersion(),x1+5,baseY +hmWidth + g.getFont().getSize());
        hmp.setLabelAbove(true);

        int track = 1;
        for (Genome genome : genomes) {
            int trackY = y2 - (int)((track + .5) * pixelsPerTrack);
            Region best = model.getBestRegion(genome);
            if (best == null) {continue;}
            String bestc = best.getChrom();
            int ostart = best.getStart();
            int oend = best.getEnd();
            hmp.setRegion(best);
            hmp.paintItem(g,x1,trackY,x2,trackY+hmWidth);
            g.setColor(Color.BLACK);
            g.drawString(genome.getVersion(),x1+5,trackY - 3);            
            java.util.List<MultiZAlignRegion> alignedRegions = model.getAlignedRegions(genome);        
            for (MultiZAlignRegion align : alignedRegions) {
                g.setColor(cs.colorAt(hashPosPair(align.getStart(),
                                                  align.getEnd(),
                                                  align.getOtherStart(),
                                                  align.getOtherEnd(),
                                                  cs.colorCount())));
                int rsx = getXPos(align.getStart(),rstart,rend,x1,x2);
                int osx = getXPos(align.getOtherStart(),ostart,oend,x1,x2);
                int rex = getXPos(align.getEnd(),rstart,rend,x1,x2);
                int oex = getXPos(align.getOtherEnd(),ostart,oend,x1,x2);

                // flip rsx and rex if on - strand
                if (align.getStrand() == '-') {
                    int temp = rsx;
                    rsx = rex;
                    rex = temp;
                }

                if (rsx >= x1 && osx >= x1) {
                    g.drawLine(rex,baseY - horizLineOffset ,oex,trackY + hmWidth + horizLineOffset);
                }
                if (rex <= x2 && oex <= x2) {
                    g.drawLine(rsx,baseY - horizLineOffset,osx,trackY + hmWidth + horizLineOffset);
                }
                g.setStroke(medium);
                g.drawLine(rsx,baseY - horizLineOffset,rex,baseY - horizLineOffset);
                g.drawLine(osx,trackY + horizLineOffset + hmWidth,
                           oex,trackY + horizLineOffset + hmWidth);
            }
            track++;
        }
        g.setStroke(oldStroke);
        g.setFont(oldFont);
        g.setColor(Color.BLACK);
    }

    private void paintMismatch(Set<Genome> genomes,
                               Genome currentGenome,
                               Graphics2D g,
                               int x1, int y1, 
                               int x2, int y2) {
        cs.reset();
        int pixelsPerTrack = (y2 -  y1) / (1 + genomes.size());
        int rstart = getRegion().getStart();
        int rend = getRegion().getEnd();
        Font oldFont = g.getFont();
        int fontSize = attrib.getFontSize(x2-x1,y2-y1);
        fontSize = Math.max(8, (int)Math.floor((double)fontSize));
        g.setFont(new Font(myFont.getFontName(), myFont.getStyle(), fontSize));
        g.setColor(Color.BLACK);
        int baseY = y2 - pixelsPerTrack/2;
        // this is the base 
        Stroke oldStroke = g.getStroke();
        Stroke thin = new BasicStroke(1);
        Stroke medium = new BasicStroke(3);
        Stroke thick = new BasicStroke(6);
        int hmWidth = 20;
        int horizLineOffset = 6;

        hmp.setLabelAbove(false);
        hmp.setRegion(getRegion());
        hmp.paintItem(g,x1,baseY,x2,baseY+hmWidth);
        g.drawString(currentGenome.getVersion(),x1+5,baseY +hmWidth + g.getFont().getSize());
        hmp.setLabelAbove(true);

        int track = 1;
        for (Genome genome : genomes) {
            int trackY = y2 - (int)((track + .5) * pixelsPerTrack);
            Region best = model.getBestRegion(genome);
            if (best == null) {continue;}
            String bestc = best.getChrom();
            int ostart = best.getStart();
            int oend = best.getEnd();
            hmp.setRegion(best);
            hmp.paintItem(g,x1,trackY,x2,trackY+hmWidth);
            g.setColor(Color.BLACK);
            g.drawString(genome.getVersion(),x1+5,trackY - 3);
            /* draw inverted chunks */
            java.util.List<MultiZAlignRegion> alignedRegions = model.getAlignedRegions(genome);        
            for (MultiZAlignRegion align : alignedRegions) {
                if (align.getStrand() != '-') {
                    continue;
                }

                g.setColor(cs.colorAt(hashPosPair(align.getStart(),
                                                  align.getEnd(),
                                                  align.getOtherStart(),
                                                  align.getOtherEnd(),
                                                  cs.colorCount())));
                int rsx = getXPos(align.getStart(),rstart,rend,x1,x2);
                int osx = getXPos(align.getOtherStart(),ostart,oend,x1,x2);
                int rex = getXPos(align.getEnd(),rstart,rend,x1,x2);
                int oex = getXPos(align.getOtherEnd(),ostart,oend,x1,x2);
                // flip rsx and rex if on - strand
                if (align.getStrand() == '-') {
                    int temp = rsx;
                    rsx = rex;
                    rex = temp;
                }
                if (rsx >= x1 && osx >= x1) {
                    g.drawLine(rex,baseY - horizLineOffset ,oex,trackY + hmWidth + horizLineOffset);
                }
                if (rex <= x2 && oex <= x2) {
                    g.drawLine(rsx,baseY - horizLineOffset,osx,trackY + hmWidth + horizLineOffset);
                }
                g.setStroke(medium);
                g.drawLine(rsx,baseY - horizLineOffset,rex,baseY - horizLineOffset);
                g.drawLine(osx,trackY + horizLineOffset + hmWidth,
                           oex,trackY + horizLineOffset + hmWidth);
                //                System.err.println("Drew inversion " + align);
            }
            /* Now draw indels */
            for (int i = 0; i < alignedRegions.size() - 1; i++) {
                MultiZAlignRegion align = alignedRegions.get(i);
                MultiZAlignRegion nextAlign = alignedRegions.get(i+1);
                g.setColor(cs.colorAt(hashPosPair(align.getStart(),
                                                  align.getEnd(),
                                                  align.getOtherStart(),
                                                  align.getOtherEnd(),
                                                  cs.colorCount())));
                int gap = nextAlign.getStart() - align.getEnd();
                int otherGap = nextAlign.getOtherStart() - align.getOtherEnd();
                if (gap != otherGap) {
                    /* if it's not a gap, it's a rearrangement */
                    boolean isgap = align.getOtherEnd() < nextAlign.getOtherEnd();
                    /* If something else falls into the apparent gap, then it's
                       not a gap but a rearrangement.
                    */
                    for (int j = 0; isgap && j < alignedRegions.size(); j++) {
                        if (j == i || j == i +1) {
                            continue;
                        }
                        MultiZAlignRegion third = alignedRegions.get(j);
                        //                        System.err.println("  Seeing if " + third + " fills the gap between " + align + " and " + nextAlign);
                        if (third.getOtherStart() >= align.getOtherEnd() &&
                            third.getOtherStart() <= nextAlign.getOtherStart()) {
                            isgap = false;
                            //                            System.err.println("    yes!");
                            break;
                        }
                    }
                    int rsx, osx, rex, oex;
                    if (isgap) {
                        rsx = getXPos(align.getEnd(),rstart,rend,x1,x2);
                        osx = getXPos(align.getOtherEnd(),ostart,oend,x1,x2);
                        rex = getXPos(nextAlign.getStart(),rstart,rend,x1,x2);
                        oex = getXPos(nextAlign.getOtherStart(),ostart,oend,x1,x2);
                        //                        System.err.println("Drew gap for " + align + " and " + nextAlign);                    
                    } else {                        
                        rsx = getXPos(align.getStart(),rstart,rend,x1,x2);
                        osx = getXPos(align.getOtherStart(),ostart,oend,x1,x2);
                        rex = getXPos(align.getEnd(),rstart,rend,x1,x2);
                        oex = getXPos(align.getOtherEnd(),ostart,oend,x1,x2);
                        //                        System.err.println("Drew rearrangement for " + align + " and " + nextAlign);                    
                    }
                    if (rsx >= x1 && osx >= x1) {
                        g.drawLine(rex,baseY - horizLineOffset ,oex,trackY + hmWidth + horizLineOffset);
                    }
                    if (rex <= x2 && oex <= x2) {
                        g.drawLine(rsx,baseY - horizLineOffset,osx,trackY + hmWidth + horizLineOffset);
                    }
                    g.setStroke(medium);
                    g.drawLine(rsx,baseY - horizLineOffset,rex,baseY - horizLineOffset);
                    g.drawLine(osx,trackY + horizLineOffset + hmWidth,
                                   oex,trackY + horizLineOffset + hmWidth);
                }

            }


            track++;
        }
        g.setStroke(oldStroke);
        g.setFont(oldFont);
        g.setColor(Color.BLACK);
    }

    public int hashPosPair (int x1, int x2, int x3, int x4, int max) {
        return Math.abs(x1 + x2 * 5 + x3 * 3 + x4 * 11) % max;
    }

}
