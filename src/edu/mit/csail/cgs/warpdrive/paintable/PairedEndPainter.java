package edu.mit.csail.cgs.warpdrive.paintable;

import java.io.File;
import java.awt.*;
import java.util.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;

import edu.mit.csail.cgs.warpdrive.model.PairedEndModel;

public class PairedEndPainter extends RegionPaintable {

    private PairedEndModel model;
    private PairedEndProperties props;
    private DynamicAttribute attrib;

    public PairedEndPainter (PairedEndModel model) {
        super();
        this.model = model;
        props = new PairedEndProperties();
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
    public PairedEndProperties getProperties() {return props;}
    public void setProperties(PairedEndProperties p) {props = p;}
    public void savePropsInDir(File dir) {
        super.savePropsInDir(dir);
        saveModelPropsInDir(dir,model);
    }
    public void loadPropsInDir(File dir) {
        super.loadPropsInDir(dir);
        loadModelPropsInDir(dir,model);
    }    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }
    public synchronized void eventRegistered(EventObject e) {        
        if (e.getSource() == model && model.isReady()) {
            setCanPaint(true);
            setWantsPaint(true);
            notifyListeners();
        }
    }
    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }
    public void paintItem(Graphics2D g, 
            int x1, int y1, 
            int x2, int y2) {
        if (!canPaint()) {
            return;
        }
        if(!model.isReady()) { return; }

        int width = x2 - x1;
        int height = Math.max(y2 - y1,1);
        int regionStart = model.getRegion().getStart();
        int regionEnd = model.getRegion().getEnd();
        int linewidth = Math.max(getProperties().LineWidth,1);
        Stroke oldStroke = g.getStroke();
        g.setStroke(new BasicStroke((float)linewidth));        
        java.util.List<PairedHit> hits = model.getResults();
        if (getProperties().DrawTrackLabel) {
            g.setFont(attrib.getLargeLabelFont(width,height));
            g.setColor(Color.BLACK);
            g.drawString("Paired " +getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
        }
        if (hits.size() == 0) { return;}
        //        int alphastep = Math.min(255, Math.max(255 / (height / (hits.size() * linewidth)), 1));
        int alphastep = 255;
        int h = height;
        int scan = 0;
        Color plusColor = new Color(0, 0, 255, alphastep);
        Color minusColor = new Color(255, 0, 0, alphastep);
        Color gray = new Color(200,200,200, alphastep);
        Color green = new Color(0, 255,0, alphastep);
        for (int i = 0; i < hits.size(); i++) {
            PairedHit hit = hits.get(i);
            int leftx1 = getXPos(hit.leftPos, regionStart, regionEnd, x1, x2);
            int leftx2 = getXPos(hit.leftStrand ? hit.leftPos + hit.leftLength : hit.leftPos - hit.leftLength,
                                 regionStart, regionEnd, x1, x2);
            if (leftx1 > leftx2) {
                int x = leftx2;
                leftx2 = leftx1;
                leftx1 = x;
            }
            int rightx1 = getXPos(hit.rightPos, regionStart, regionEnd, x1, x2);
            int rightx2 = getXPos(hit.rightStrand ? hit.rightPos + hit.rightLength : hit.rightPos - hit.rightLength,
                                  regionStart, regionEnd, x1, x2);            
            if (rightx1 > rightx2) {
                int x = rightx2;
                rightx2 = rightx1;
                rightx1 = x;
            }
            g.setColor(hit.leftStrand == hit.rightStrand ? green : gray);
            if (leftx2 < rightx1) {
                g.drawLine(leftx2, y1+h, rightx1, y1+h);
            } else {
                g.drawLine(rightx2, y1+h, leftx1, y1+h);
            }
            if (leftx2 == leftx1) {leftx2++;}
            g.setColor(hit.leftStrand ? plusColor : minusColor);
            g.drawLine(leftx1, y1 + h, leftx2, y1+h);
            if (rightx2 == rightx1) {rightx2++;}
            g.setColor(hit.rightStrand ? plusColor : minusColor);
            g.drawLine(rightx1, y1+h, rightx2, y1+h);

            
            h -= linewidth * 2;
            if (h < 0) {
                h = height;
            }            
        }
        g.setStroke(oldStroke);
    }

}
