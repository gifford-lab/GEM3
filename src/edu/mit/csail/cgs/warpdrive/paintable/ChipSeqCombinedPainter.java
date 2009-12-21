package edu.mit.csail.cgs.warpdrive.paintable;

import edu.mit.csail.cgs.warpdrive.model.ChipSeqDataModel;
import java.awt.Graphics2D;

public class ChipSeqCombinedPainter extends ChipSeqPainter {

    private ChipSeqBasicOverlapPainter basic;
    private ChipSeqAboveBelowStrandPainter stranded;

    public ChipSeqCombinedPainter (ChipSeqDataModel model) {
        super(model);
        basic = new ChipSeqBasicOverlapPainter(model);
        basic.setProperties(getProperties());
        stranded = new ChipSeqAboveBelowStrandPainter(model);
        stranded.setProperties(getProperties());       
    }

    protected void paintOverlapping(Graphics2D g, 
                                    int x1, int y1, 
                                    int x2, int y2) {
        if (getProperties().Stranded) {
            stranded.paintOverlapping(g,x1,y1,x2,y2);
        } else {
            basic.paintOverlapping(g,x1,y1,x2,y2);
        }
    }
    
    protected void paintNonOverlapping(Graphics2D g, 
                                       int x1, int y1, 
                                       int x2, int y2) {
        if (getProperties().Stranded) {
            stranded.paintNonOverlapping(g,x1,y1,x2,y2);
        } else {
            basic.paintNonOverlapping(g,x1,y1,x2,y2);
        }
    }
}