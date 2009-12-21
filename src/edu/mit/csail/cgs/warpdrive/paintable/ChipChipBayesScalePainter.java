package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.LineMetrics;
import java.util.*;
import java.awt.geom.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;

public class ChipChipBayesScalePainter extends ChipChipScalePainter {

    static String problabel = "Probability", strlabel = "Strength", jbdlabel = "Binding Prob";
    private ChipChipBayesProperties props;

    public ChipChipBayesScalePainter (ChipChipScaleModel m, Listener<EventObject> parent, ChipChipBayesProperties p) {
        super(m,parent,p);
        props = p;
    }
    public ChipChipBayesProperties getProperties() {return props;}
    public void setProperties(ChipChipBayesProperties p) {props = p;}

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        Font oldFont = g.getFont();
        g.setFont(DynamicAttribute.getGlobalAttributes().getRegionLabelFont(lrx-ulx,lry - uly));
        FontMetrics fontmetrics = g.getFontMetrics();
        g.setColor(Color.BLACK);
        g.drawLine(ulx,lry,lrx,lry);

            AffineTransform oldtransform = g.getTransform();
            LineMetrics linemetrics = fontmetrics.getLineMetrics(problabel,g);
            float sh = linemetrics.getHeight();
            float sw = fontmetrics.charsWidth(strlabel.toCharArray(),0,strlabel.length());
            linemetrics = fontmetrics.getLineMetrics(problabel,g);
            float ph = linemetrics.getHeight();
            float pw = fontmetrics.charsWidth(problabel.toCharArray(),0,problabel.length());

            super.paintItem(g,ulx,uly,lrx - (int)(ph * 1.2),lry);

            g.rotate(Math.PI / -2,0,0);        
            g.setColor(Color.WHITE);
            g.fillRect(-1 * (int)((lry + uly) / 2 + sw / 1.9), (int)(sh * 1.1), (int)(sw*1.2),(int)(sh * 1.2));
            g.setColor(Color.GRAY);
            g.drawString(strlabel, -1 * (int)((lry + uly) / 2 + sw / 2), (int)(sh * 1.1));
            g.setColor(Color.BLUE);
            g.drawString(problabel,-1 * (int)((lry + uly) / 2 + pw / 2), lrx - 2);
            g.setTransform(oldtransform);
            int y = getYPos(.25,0,1,uly,lry,false) + (int)(ph/2);
            g.drawString("0.25",lrx - fontmetrics.charsWidth("0.25".toCharArray(),0,4),y);
            y = getYPos(.75,0,1,uly,lry,false) + (int)(ph/2);
            g.drawString("0.75",lrx - fontmetrics.charsWidth("0.75".toCharArray(),0,4),y);
            y = getYPos(1.0,0,1,uly,lry,false) + (int)(ph/2);
            g.drawString("1.0",lrx - fontmetrics.charsWidth("1.0".toCharArray(),0,3),y);
           
        g.setFont(oldFont);

    }
}
