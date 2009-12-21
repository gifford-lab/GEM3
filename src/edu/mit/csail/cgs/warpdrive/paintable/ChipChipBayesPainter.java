package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;

public class ChipChipBayesPainter extends ChipChipPainter {

    private ChipChipBayes data;
    private ChipChipBayesProperties props;

    public ChipChipBayesPainter (ChipChipBayes data, ChipChipDataModel model) {
        super(data,model);
        this.data = data;
        props = new ChipChipBayesProperties();
    }

    public ChipChipBayesProperties getProperties() {
        return props;
    }
    public void setProperties(ChipChipBayesProperties p) {props = p;}

    public void paintItem(Graphics2D g, 
                          int ulx, int uly, 
                          int lrx, int lry) {
        if (!canPaint()) {
            return;
        }
        int i = 0;
        // get color from track properties or pick one
        int startindex = data.baseToIndex(getRegion().getStart());
        int endindex = data.baseToIndex(getRegion().getEnd());
        //        System.err.println("startindex is " + startindex + " and end is " + endindex);

        int lastx = -1;
        int lasty = -1;
        int lastpy = -1;
        for (i = startindex; i < endindex; i++) {
            double maxratio;
            if (getProperties().MaxRatio > 0) {
                maxratio = getProperties().MaxRatio;
            } else {
                maxratio = getScaleModel().getMaxVal();
            }
            int x = getXPos(data.getPos(i),
                            getRegion().getStart(), getRegion().getEnd(),
                            ulx,lrx);
            int y,py;
            y = getYPos(data.getValue(i,1),
                        0,
                        maxratio,
                        uly,lry,
                        getProperties().RatiosOnLogScale);
            py = getYPos(data.getPosterior(i),
                         0.0,
                         1.0,
                         uly,lry,
                         getProperties().RatiosOnLogScale);                
            g.setColor(Color.GRAY);
            paintDatapointAt(g,x,y,i,1);            
            if ((lastx != -1) && 
                (Math.abs(data.getPos(i) - data.getPos(i-1)) < 500)) {
                connectDatapoints(g,lastx,lasty,x,y);
            }
            g.setColor(Color.BLUE);
            paintDatapointAt(g,x,py,i,1);
            if ((lastx != -1) && 
                (Math.abs(data.getPos(i) - data.getPos(i-1)) < 500)) {
                connectDatapoints(g,lastx,lastpy,x,py);
            }
            lastx = x;
            lasty = y;
            lastpy = py;
        }
        if (props.DrawTrackLabel) {
            Font oldFont = g.getFont();
            g.setFont(DynamicAttribute.getGlobalAttributes().getLargeLabelFont(lrx-ulx,lry-uly));
            g.setColor(Color.BLACK);
            g.drawString(getLabel(), ulx + 5, uly + g.getFont().getSize() + 5);
            g.setFont(oldFont);
        }
    }
}

