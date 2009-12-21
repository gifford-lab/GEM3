/*
 * Created on Mar 28, 2006
 */
package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import java.text.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.nouns.*;

/**
 * @author tdanford
 */
public class FixedWidthHashPaintable extends RegionPaintable {
    
    private static NumberFormat nf;
    private boolean labelAbove;
    private PaintableProperties props;
    static {
        nf = NumberFormat.getInstance();
    }    

    public FixedWidthHashPaintable () {
        super();
        labelAbove = true;
        props = new PaintableProperties();
    }
    public PaintableProperties getProperties() {return props;}

    public void setLabelAbove (boolean b) {labelAbove = b;}

    public boolean canPaint() {
        return true;
    }

    public void setRegion(Region r) {
        super.setRegion(r);
        setCanPaint(true);
        setWantsPaint(true);
        notifyListeners();
    }

    public int getMaxVertSpace() { return 30; }
    
    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.viz.paintable.Paintable#paintItem(java.awt.Graphics, int, int, int, int)
     */
    public void paintItem(Graphics2D g, int x1, int y1, int x2, int y2) {
    	Region region = getRegion();
        int rs = getRegion().getStart(), re = getRegion().getEnd();
        int rw = re - rs + 1;
        int h = y2 - y1;
        int w = x2 - x1;
        
        double log10 = Math.log(10.0);
        double logScale = Math.log((double)rw) / log10;
        Font oldFont = g.getFont();
        g.setFont(new Font("Arial",Font.PLAIN,(int)(h * .4)));
        FontMetrics fontmetrics = g.getFontMetrics();
        char[] nines = {'9','9','9','9','9','9','9','9','9','9','9','9','9'};


        //int tickSize = Math.max(1, pow10((int)Math.floor(logScale) - 2));
        
        int tickSize = 10000;
        int tickDigits = (int)Math.ceil(Math.log(re / tickSize) / log10) + 1;
        int spaceUsed = fontmetrics.charsWidth(nines,0,tickDigits);
        
        double multiplier = 1.2;
        int spaceAvail = (int) (multiplier * w / (rw / ((double)tickSize)));

        /*
        while(spaceUsed > spaceAvail) {
            tickSize *= 10;
            tickDigits = (int)Math.ceil(Math.log(re / tickSize) / log10) + 1;
            spaceUsed = fontmetrics.charsWidth(nines,0,tickDigits);
            spaceAvail = (int) (multiplier * w / (rw / ((double)tickSize)));
        }

        tickSize = Math.max(tickSize, 1);
        */
        
        int hh = h/2;
        int my, ly;
        if (labelAbove) {
            my = y1 + hh;
            ly = my - 2;
        } else {
            my = y1;
            ly = y2;
        }
        
        Stroke oldStroke = g.getStroke();
        g.setStroke(new BasicStroke((float)2.0));

        g.setColor(Color.white);
        g.fillRect(x1, my, w, hh);
        
        g.setColor(Color.black);
        //g.drawLine(x1, my, x2, my);
        
        int bpStart = tickSize * (int)Math.ceil((double)rs / (double)tickSize);
        double scale = (double)w / (double)rw;
        int tickPixWidth = (int)Math.round((double)tickSize * scale);
        
        boolean paintBlack = true;
        Region r = getRegion();
        String axisLabel = String.format("%s: %s-%s,   %s bp/segment", 
                r.getChrom(), nf.format(r.getStart()), nf.format(r.getEnd()), 
                nf.format(tickSize));
        if (!axisLabel.matches("^chr.*")) {
            axisLabel = "chr" + axisLabel;
        }
        
        /*
        String axisLabel = String.format("%s bp/segment", 
                nf.format(tickSize));
         */
        
        g.setColor(Color.black);
        //g.drawString(axisLabel, x1 + 2, ly);

        int leftbp = rs - rs%tickSize + tickSize;
        int rightbp = re - re%tickSize;
        
        String leftstring = String.format("%d", leftbp);
        String rightstring = String.format("%d", rightbp);
        int rightwidth = fontmetrics.stringWidth(rightstring);
        
        g.drawString(leftstring, x1+2, my-2);
        g.drawString(rightstring, x2-rightwidth, my-2);

        for(int bp = bpStart; bp <= re; bp += tickSize) { 
            int xOffset = x1 + (int)Math.round((double)(bp-rs) * scale);
        
            g.setColor(Color.black);
            //g.drawLine(xOffset, my, xOffset, my+hh);   
            g.drawLine(xOffset, y1, xOffset, my);   
            if(paintBlack) { 
            	//g.drawString(String.valueOf(bp / tickSize), xOffset+2, my + hh - 1);
            }

            paintBlack = !paintBlack;
        }
        
        g.setFont(oldFont);
        g.setStroke(oldStroke);
    }
    
    private int pow10(int exp) { 
        int value = 1;
        for(int i = 0; i < exp; i++) { 
            value *= 10;
        }
        return value;
    }
}
