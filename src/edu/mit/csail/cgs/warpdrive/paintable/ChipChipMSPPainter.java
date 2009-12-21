package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMSP;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;

public class ChipChipMSPPainter extends TimChipChipPainter {
    private ChipChipMSP data;

    public ChipChipMSPPainter (ChipChipMSP data,  ChipChipDataModel model) {
        super(data,model);
        this.data = data;
    }    

    public void paintDataPointAt(Graphics2D g, int x,int y,int i,int j) {
        
        Color pc = g.getColor();
        
        int squaresize = 4;
        if (data.getPval3(i) < .001) {
            squaresize = 8;
        } else if (data.getPval3(i) < .01) {
            squaresize = 6;
        }
        
        //g.setColor(Color.white);
        //g.fillRect(x-squaresize/2,y-squaresize/2,squaresize,squaresize);
        
        //g.setColor(pc);
        //g.drawRect(x-squaresize/2,y-squaresize/2,squaresize,squaresize);
    }
}

