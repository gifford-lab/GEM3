package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseEvent;
import java.sql.SQLException;
import java.util.*;
import javax.swing.JToolTip;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.warpdrive.components.GOAnnotationPanel;
import edu.mit.csail.cgs.warpdrive.model.RegionExpanderModel;

public class HarbisonProbePainter extends RegionPaintable {
    
    private RegionExpanderModel<HarbisonRegCodeProbe> model;
    private Font myFont;
    private DynamicAttribute attrib;
    private Vector<HarbisonRegCodeProbe> probes;    
    private PaintableProperties props;
    
    public HarbisonProbePainter(RegionExpanderModel<HarbisonRegCodeProbe> model) {
        super();
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        myFont = new Font("Arial", Font.PLAIN, 12);
        model.addEventListener(this);
        props = new PaintableProperties();
        initLabels();
    }
    public PaintableProperties getProperties() {return props;}
    
    public void clickedOnItem(ActionEvent e) {
        System.err.println(e);
        String iName = e.getActionCommand();        
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public synchronized void eventRegistered(EventObject e) {        
        if (e.getSource() == model &&
            model.isReady()) {            

            setCanPaint(true);
            setWantsPaint(true);
            probes = null;
            setLayoutProbes();
            notifyListeners();
        }
    }    

    public int getMaxVertSpace() { 
        return 25;
    }

    private void setLayoutProbes() {
        if (canPaint() && probes == null) {
            Iterator<HarbisonRegCodeProbe> itr = model.getResults();
            probes = new Vector<HarbisonRegCodeProbe>();
            while(itr.hasNext()) { probes.add(itr.next()); }
        }
    }
    
    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        if (!canPaint()) {
            return;
        }
               
        int w = x2 - x1, h = y2 - y1;
        int my = y1 + (h / 2);
 
        Region region = model.getRegion();
        int rs = region.getStart(), re = region.getEnd();
        int rw = re - rs + 1;

        // Additions by Lindy
        // Dynamically change font size as the size of the window changes
        int fontSize = attrib.getFontSize(w,h);
        fontSize = Math.max(8, (int)Math.floor((double)fontSize * 0.90));
        Font oldFont = g.getFont();
        g.setFont(new Font(myFont.getFontName(), myFont.getStyle(), fontSize));
        //--------------------------------------------------
        
        double xScale = (double)w / (double)rw;
        clearLabels();
        int fontsize = g.getFont().getSize();

        for(HarbisonRegCodeProbe probe : probes) { 
            int ps = probe.getStart(), pe = probe.getEnd();
            int px1 = x1 + (int)Math.round(xScale * (double)(ps-rs));
            int px2 = x1 + (int)Math.round(xScale * (double)(pe-rs));
            int py1 = y1;
            int pw = px2 - px1;
            int ph = h;
            
            g.setColor(Color.gray);
            if(probe.size() > 0) { 
                g.fillRect(px1, py1, pw, ph);                                
            } else { 
                g.drawRect(px1, py1, pw, ph);                
            }
            addLabel(px1, py1, pw, ph, probe.getName());
            
            for(String f : probe.getFactors()) { 
                for(String c : probe.getConditions(f)) { 
                    String lbl = f + "_" + c + " (" + 
                    HarbisonRegCodeProbe.strengthStrings[probe.getBindingStrength(f, c)] + ")";
                    addLabel(px1, py1, pw, ph, lbl);
                }
            }
            
            if (probe.getName().length() * fontsize < pw && fontsize < ph) {
                g.setColor(Color.black);
                g.drawString(probe.getName(), px1 + 2, py1 + ph - 2);
            }
        }

        g.setFont(oldFont);
        if (props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1,y2);
        }
        
    }

}
