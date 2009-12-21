package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.GenericExperiment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.ewok.nouns.*;

public class ChipChipIntensitiesPainter extends TimChipChipPainter {
    private ChipChipDataModel model;
    public ChipChipIntensitiesPainter (ChipChipData data, ChipChipDataModel model) {
        super(data,model);
        this.model = model;
    }
    
    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        //        super.paintItem(g,x1,y1,x2,y2);        
        if (model.getGenericExperiment() instanceof ChipChipData) {
            double maxval = 65535;
            Region region = model.getRegion();
            int rs = region.getStart(), re = region.getEnd();
            ChipChipData data = (ChipChipData)model.getGenericExperiment();
            int circleradius = 1;
            int circlediameter = 1;
            for (int i = 0; i < data.getCount(); i++) {
                int px = getXPos(data.getPos(i),
                                 rs,re,x1,x2);
                for (int j = 0; j < data.getReplicates(i); j++) {
                    int ippy = getYPos(data.getIP(i,j),
                                       0,maxval,
                                       y1,y2,getProperties().IntensitiesOnLogScale);
                    int wcepy = getYPos(data.getWCE(i,j),
                                        0,maxval,
                                        y1,y2,getProperties().IntensitiesOnLogScale);
                    g.setColor(Color.GREEN);
                    g.fillOval(px-circleradius,ippy-circleradius,circlediameter,circlediameter);
                    g.setColor(Color.RED);
                    g.fillOval(px-circleradius,wcepy-circleradius,circlediameter,circlediameter);
                }
            }
        }
        if (getProperties().DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1,y1 + g.getFont().getSize() * 2);    	
        }
    }
}

