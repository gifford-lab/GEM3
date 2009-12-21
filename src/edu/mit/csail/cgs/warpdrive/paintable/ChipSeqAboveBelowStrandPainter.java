package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.font.LineMetrics;
import java.awt.geom.Rectangle2D;
import java.util.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.GenericExperiment;
import edu.mit.csail.cgs.datasets.chippet.WeightedRunningOverlapSum;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.ExonicGene;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.ChipSeqDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipSeqScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.ewok.nouns.*;

public class ChipSeqAboveBelowStrandPainter extends ChipSeqPainter {

    private Color color;
    protected Vector<ChipSeqHit> watsonLayoutHits = new Vector<ChipSeqHit>();
    protected Vector<ChipSeqHit> crickLayoutHits = new Vector<ChipSeqHit>();
    protected NonOverlappingLayout<ChipSeqHit> watsonLayout = new NonOverlappingLayout<ChipSeqHit>();
    protected NonOverlappingLayout<ChipSeqHit> crickLayout = new NonOverlappingLayout<ChipSeqHit>();


    public ChipSeqAboveBelowStrandPainter(ChipSeqDataModel model) {
        super(model);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
    
    protected void setLayoutHits() {        
        if (canPaint()) {        	
            if (getProperties().Overlapping) {
                /* don't do anything.  The model can provide us with the RunningOverlapSums */
            } else {
                Iterator<ChipSeqHit> itr = model.getResults();
                watsonLayoutHits = new Vector<ChipSeqHit>();
                crickLayoutHits = new Vector<ChipSeqHit>();
                Region extended = null;
                while (itr.hasNext()) {
                    ChipSeqHit hit = itr.next();
                    if (hit.getStrand() == '+') {
                        watsonLayoutHits.add(hit);
                    } else {
                        crickLayoutHits.add(hit);
                    }
                }
                watsonLayout.setRegions(watsonLayoutHits);
                crickLayout.setRegions(crickLayoutHits);
            }
        }
    }
    
    
    protected void paintOverlapping(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        int width = x2 - x1, height = y2 - y1;
        Region r = model.getRegion();
                
        WeightedRunningOverlapSum watsonOverlap = model.getWatsonRunningOverlap();
        WeightedRunningOverlapSum crickOverlap = model.getCrickRunningOverlap();
        double watsonMaxOverlap = watsonOverlap.getMaxOverlap();
        double crickMaxOverlap = crickOverlap.getMaxOverlap();
        
        int baseline = y1 + height/2;
        int availHeight = (int)Math.round(0.475 * height);
        
        g.setColor(Color.black);
        g.drawLine(x1, baseline, x2, baseline);

        double maxobsoverlap = Math.max(watsonMaxOverlap, crickMaxOverlap);
        int propmaxoverlap = ((ChipSeqProperties)getProperties()).MaxReadCount;
        /* maxoverlap is the maximum height of the reads scale
           If the property is negative, use the observed maximum (ie, dynamic scaling).
           Otherwise, use the value from the property and clip peaks
           higher than that value
        */
        double maxoverlap;
        if (propmaxoverlap > 0) {
            maxoverlap = propmaxoverlap;
        } else {
            maxoverlap = maxobsoverlap;
        }
        maxoverlap = Math.max(maxoverlap, 4);
        
        double trackHeight = (double)availHeight / (double)(maxoverlap+1.0);
        int pixTrackHeight = Math.max(1, (int)Math.floor(trackHeight));
        
        g.setColor(Color.blue);
        Stroke oldStroke = g.getStroke();
        //g.setStroke(new BasicStroke((float)2.0));

        for(int over = 1; over <= Math.min(Math.ceil(watsonMaxOverlap), Math.ceil(maxoverlap)); over++) { 
            Collection<Region> watsonRegs = watsonOverlap.collectRegions(over);
            for(Region or : watsonRegs) { 
                int orx1 = calcX(or.getStart(), r.getStart(), r.getEnd(), x1, x2);
                int orx2 = calcX(or.getEnd(), r.getStart(), r.getEnd(), x1, x2);
                int ory1 = baseline - 1 - (int)Math.floor(trackHeight * (double)over);
                g.drawRect(orx1, ory1, (orx2-orx1), pixTrackHeight);
            }
        }

        g.setColor(Color.red);
        for(int over = 1; over <= Math.min(Math.ceil(crickMaxOverlap), Math.ceil(maxoverlap)); over++) { 
            Collection<Region> crickRegs = crickOverlap.collectRegions(over);
            for(Region or : crickRegs) { 
                int orx1 = calcX(or.getStart(), r.getStart(), r.getEnd(), x1, x2);
                int orx2 = calcX(or.getEnd(), r.getStart(), r.getEnd(), x1, x2);
                int ory1 = baseline + 1 - pixTrackHeight + (int)Math.floor(trackHeight * (double)over);
                g.drawRect(orx1, ory1, (orx2-orx1), pixTrackHeight);
            }
        }
        /* draw the scale */
        g.setColor(Color.black);
        g.setFont(attrib.getLargeLabelFont(width,height));
        int step = Math.max(1,(int)Math.round(maxoverlap / 5));
        for (int i = step; i <= Math.ceil(maxoverlap); i += step) {
            int ypos = baseline + 1 - pixTrackHeight + (int)Math.floor(trackHeight * (double)i);
            g.drawString(Integer.toString(i),
                         5,
                         ypos);
            ypos = baseline - 1 - (int)Math.floor(trackHeight * (double)i);
            g.drawString(Integer.toString(i),
                         5,
                         ypos);
        }        
        
        g.setStroke(oldStroke);
    }
    
    
    protected void paintNonOverlapping(Graphics2D g, 
            int x1, int y1, 
            int x2, int y2) {

        int width = x2 - x1;
        int height = y2 - y1;
        
        int baseline = y1 + height/2;
        int availHeight = (int)Math.round(0.475 * height);
        
        g.setColor(Color.black);
        g.drawLine(x1, baseline, x2, baseline);

        /* if the user has asked for a maxcount greater than
           the maximum number of tracks into which NonOverlappingLayout
           would have done the layout, tell NOL to use more tracks
        */
        int propmaxTracks = ((ChipSeqProperties)getProperties()).MaxReadCount;
        if (watsonLayout.getMaxTracks() <= propmaxTracks) {
            watsonLayout.setMaxTracks(propmaxTracks);
        }
        if (crickLayout.getMaxTracks() <= propmaxTracks) {
            crickLayout.setMaxTracks(propmaxTracks);
        }
        int watsonNumTracks = watsonLayout.getNumTracks();
        int crickNumTracks = crickLayout.getNumTracks();
        int maxobsTracks = Math.max(watsonNumTracks, crickNumTracks);
        /* maxoverlap is the maximum height of the reads scale
           If the property is negative, use the observed maximum (ie, dynamic scaling).
           Otherwise, use the value from the property and clip peaks
           higher than that value
        */
        int numTracks;
        if (propmaxTracks > 0) {
            numTracks = propmaxTracks;
        } else {
            numTracks = maxobsTracks;
        }
        numTracks = Math.max(numTracks, 4);
        
        double trackHeight = (double)availHeight / (double)(numTracks+1.0);
        int pixTrackHeight = Math.max(2, (int)Math.floor(trackHeight));
        int halfTrackHeight = Math.max(1, pixTrackHeight/2);
        
        Region region = model.getRegion();
        int rs = region.getStart();
        int re = region.getEnd();
        int rw = re - rs + 1;
        double xScale = (double)width / (double)rw;

        int hitHeight = Math.max(2, (int)Math.floor((double)pixTrackHeight * 0.90));
        int halfHitHeight = hitHeight / 2;
        
        g.setColor(Color.blue);
        Stroke oldStroke = g.getStroke();
        //g.setStroke(new BasicStroke((float)2.0));

        for(ChipSeqHit watsonHit : watsonLayoutHits) {
            int track = 0;
            
            if(!watsonLayout.hasTrack(watsonHit)) { 
                System.err.println("No track assigned to hit: " + watsonHit.getLocationString());
            } 
            else { 
                track = Math.min(numTracks,watsonLayout.getTrack(watsonHit));
            }

            int hy1 = baseline - (pixTrackHeight * (track + 1));
            int hy2 = hy1 + pixTrackHeight;
            int hmy = hy1 + halfTrackHeight;
                        
            int htop = hmy - halfHitHeight;
            int hbottom = hmy + halfHitHeight;
            
            int hitStart = watsonHit.getStart();
            int hitEnd = watsonHit.getEnd();
            
            int hx1 = xcoord(hitStart, x1, rs, xScale);
            int hx2 = xcoord(hitEnd, x1, rs, xScale);
            int hleft = Math.max(x1, hx1);
            int hright = Math.min(x2, hx2);

            int rectwidth = hright - hleft + 1;

            if (watsonHit.getWeight() == 1) {
              g.setColor(Color.blue);
            }
            else {              
              g.setColor(Color.cyan);
            }
            g.drawRect(hleft, htop, rectwidth, hbottom - htop);
        }
        
        
        g.setColor(Color.red);        
        for(ChipSeqHit crickHit : crickLayoutHits) {
            int track = 0;
            
            if(!crickLayout.hasTrack(crickHit)) { 
                System.err.println("No track assigned to hit: " + crickHit.getLocationString());
            } 
            else { 
                track = Math.min(numTracks,crickLayout.getTrack(crickHit));
            }

            int hy1 = baseline + (pixTrackHeight * track);
            int hy2 = hy1 + pixTrackHeight;
            int hmy = hy1 + halfTrackHeight;           
            
            int htop = hmy - halfHitHeight;
            int hbottom = hmy + halfHitHeight;
            
            int hitStart = crickHit.getStart();
            int hitEnd = crickHit.getEnd();
            
            int hx1 = xcoord(hitStart, x1, rs, xScale);
            int hx2 = xcoord(hitEnd, x1, rs, xScale);
            int hleft = Math.max(x1, hx1);
            int hright = Math.min(x2, hx2);

            int rectwidth = hright - hleft + 1;

            if (crickHit.getWeight() == 1) {
              g.setColor(Color.red);
            }
            else {              
              g.setColor(Color.pink);
            }
            g.drawRect(hleft, htop, rectwidth, hbottom - htop);
        }
        
        /* draw the scale */
        g.setColor(Color.black);
        g.setFont(attrib.getLargeLabelFont(width,height));
        int step = Math.max(1,numTracks / 5);
        for (int i = step; i <= numTracks; i += step) {
            int ypos = baseline - (pixTrackHeight * (i + 1));
            g.drawString(Integer.toString(i),
                         5,
                         ypos);
            ypos = baseline + (pixTrackHeight * i);
            g.drawString(Integer.toString(i),
                         5,
                         ypos);
        }

        g.setStroke(oldStroke);
    }
}
