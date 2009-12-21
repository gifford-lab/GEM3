/*
 * Created on May 11, 2006
 */
package edu.mit.csail.cgs.metagenes;

import java.io.*;
import java.util.*;
import java.awt.*;

import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.paintable.*;

/**
 * @author tdanford
 */
public class FauxGenePaintable extends AbstractPaintable {
    
    private double widthFraction;
    private DynamicAttribute attrib;

    public FauxGenePaintable(double wf) {
        widthFraction = wf;
        attrib = new DynamicAttribute();
    }

    public FauxGenePaintable(int h, int w, double wf) {
        super(h, w);
        widthFraction = wf;
        attrib = new DynamicAttribute();
    }
    
    public void setWidthFraction(double wf) { widthFraction = wf; } 

    public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
    	
        Graphics2D g2 = (Graphics2D)g;
        int w = x2 - x1, h = y2 - y1;
        int mx = x1 + (int)Math.round((double)w * widthFraction);
        
        //int lineWidth = Math.max(1, 2 * attrib.getLineWidth(w, h) / 3);
        int lineWidth = 3;
        Stroke oldStroke = g2.getStroke();
        g2.setStroke(new BasicStroke((float)lineWidth));

        int spacing = h/6 + (lineWidth * 2);
        int ah = h - spacing * 2;
        int my = y1 + h/2;
        int ty = my - ah/2;
        int by = my + ah/2;
                
        g2.setColor(Color.black);
        g2.drawLine(x1, my, mx, my);
        g2.drawLine(mx, ty, mx, by);
        g2.drawLine(mx, ty, x2, ty);
        g2.drawLine(mx, by, x2, by);

        int lxa = mx - spacing;
        int rxa = mx + spacing * 2;
        int txa = rxa + spacing;
        int[] ax = {lxa, lxa, rxa, rxa, txa, rxa, rxa }; 
        
        int bya = my;
        int tya = ty - spacing;
        int[] ay = { bya, tya, tya, tya - spacing, tya, tya + spacing, tya };
        
        g2.drawPolyline(ax, ay, 7);
        
        g2.setStroke(oldStroke);
    }
}
