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

public class ChipSeqBasicOverlapPainter extends ChipSeqPainter {

    private Color color;

    public ChipSeqBasicOverlapPainter(ChipSeqDataModel model) {
        super(model);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
    
    protected void paintOverlapping(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        int width = x2 - x1, height = y2 - y1;
        Region r = model.getRegion();
        
        WeightedRunningOverlapSum overlap = model.getTotalRunningOverlap();
        int maxTracksOverlap = (int)Math.ceil(overlap.getMaxOverlap());

        int propmaxoverlap = ((ChipSeqProperties)getProperties()).MaxReadCount;
        /* maxoverlap is the maximum height of the reads scale
           If the property is negative, use the observed maximum (ie, dynamic scaling).
           Otherwise, use the value from the property and clip peaks
           higher than that value
        */
        int maxOverlap;
        if (propmaxoverlap > 0) {
            maxOverlap = propmaxoverlap;
        } else {
            maxOverlap = maxTracksOverlap;
        }
        maxOverlap = Math.max(maxOverlap, 4);

        int baseline = y2 - (y2-y1)/5;
        int availHeight = baseline - y1;
        
        g.setColor(Color.black);
        g.drawLine(x1, baseline, x2, baseline);
        double trackHeight = (double)availHeight / (double)(maxOverlap+1);
        int pixTrackHeight = Math.max(1, (int)Math.floor(trackHeight));
        
        g.setColor(Color.blue);
        Stroke oldStroke = g.getStroke();
        //g.setStroke(new BasicStroke((float)2.0));
        
        for(int over = 1; over <= Math.min(maxOverlap, maxTracksOverlap); over++) { 
            Collection<Region> regs = overlap.collectRegions(over);
            for(Region or : regs) { 
                int orx1 = calcX(or.getStart(), r.getStart(), r.getEnd(), x1, x2);
                int orx2 = calcX(or.getEnd(), r.getStart(), r.getEnd(), x1, x2);
                int ory1 = baseline - (int)Math.floor(trackHeight * (double)over);
                g.drawRect(orx1, ory1, (orx2-orx1), pixTrackHeight);
            }
        }
        
        System.out.println(String.format("%s --> %f overlap", r.getLocationString(), model.getTotalMaxOverlap()));
        
        /* draw the scale */
        g.setColor(Color.black);
        g.setFont(attrib.getLargeLabelFont(width,height));
        int step = Math.max(1,maxOverlap / 5);
        for (int i = step; i <= maxOverlap; i += step) {
            int ypos = baseline - (int)Math.floor(trackHeight * (double)i);
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
        
        int baseline = y2 - (y2-y1)/5;
        int availHeight = baseline - y1;
        
        g.setColor(Color.black);
        g.drawLine(x1, baseline, x2, baseline);

        int numTracks = Math.max(1, totalLayout.getNumTracks());
        int propmaxoverlap = ((ChipSeqProperties)getProperties()).MaxReadCount;
        /* maxoverlap is the maximum height of the reads scale
           If the property is negative, use the observed maximum (ie, dynamic scaling).
           Otherwise, use the value from the property and clip peaks
           higher than that value
        */
        if (propmaxoverlap > 0) {
            numTracks = propmaxoverlap;
        } 
        numTracks = Math.max(numTracks, 4);

        double trackHeight = (double)availHeight / (double)(numTracks+1);
        int pixTrackHeight = Math.max(2, (int)Math.floor(trackHeight));
        int halfTrackHeight = Math.max(1, pixTrackHeight/2);
        
        Region region = model.getRegion();
        int rs = region.getStart(), re = region.getEnd();
        int rw = re - rs + 1;
        double xScale = (double)width / (double)rw;

        g.setColor(Color.blue);
        Stroke oldStroke = g.getStroke();
        //g.setStroke(new BasicStroke((float)2.0));
                        
        for(Region hit : totalLayoutHits) {
            int track = 0;
            
            if(!totalLayout.hasTrack(hit)) { 
                System.err.println("No track assigned to hit: " + hit.getLocationString());
            } else { 
                track = Math.min(numTracks, totalLayout.getTrack(hit));
            }

            int hy1 = baseline - (pixTrackHeight * (track + 1));
            int hy2 = hy1 + pixTrackHeight;
            int hmy = hy1 + halfTrackHeight;
            
            int hitHeight = Math.max(2, (int)Math.floor((double)pixTrackHeight * 0.90));
            int halfHitHeight = hitHeight / 2;
            
            int htop = hmy - halfHitHeight;
            int hbottom = hmy + halfHitHeight;
            
            int hitStart = hit.getStart();
            int hitEnd = hit.getEnd();
            
            int hx1 = xcoord(hitStart, x1, rs, xScale);
            int hx2 = xcoord(hitEnd, x1, rs, xScale);
            int hleft = Math.max(x1, hx1);
            int hright = Math.min(x2, hx2);

            int rectwidth = hright - hleft + 1;

            g.drawRect(hleft, htop, rectwidth, hbottom - htop);
        }
        g.setStroke(oldStroke);
    }

}
