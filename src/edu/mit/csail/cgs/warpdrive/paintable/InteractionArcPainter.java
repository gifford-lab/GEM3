package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.io.File;
import java.util.EventObject;

import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.warpdrive.model.InteractionArcModel;
import edu.mit.csail.cgs.warpdrive.model.PairedEndModel;

public class InteractionArcPainter extends RegionPaintable {
	
	private InteractionArcModel model;
	private InteractionArcProperties props;
	private DynamicAttribute attrib;
	
	public InteractionArcPainter (InteractionArcModel model) {
        super();
        this.model = model;
        props = new InteractionArcProperties();
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
	
	public InteractionArcProperties getProperties() {return props;}
    public void setProperties(InteractionArcProperties p) {props = p;}
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

	public void paintItem(Graphics2D g, int x1, int y1, int x2, int y2) {
		if (!canPaint()) {
            return;
        }
        if(!model.isReady()) { return; }
        
        int width = x2 - x1;
        int height = Math.max(y2 - y1,1);
        int regionStart = model.getRegion().getStart();
        int regionEnd = model.getRegion().getEnd();
        int regionWidth = model.getRegion().getWidth();
        int linewidth = Math.max(getProperties().LineWidth,1);
        Stroke oldStroke = g.getStroke();
        g.setStroke(new BasicStroke((float)linewidth));        
        java.util.List<PairedHit> hits = model.getResults();
        if (getProperties().DrawTrackLabel) {
            g.setFont(attrib.getLargeLabelFont(width,height));
            g.setColor(Color.BLACK);
            g.drawString("Paired " +getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
        }
        int h = height;
        for (int i = 0; i < hits.size(); i++) {
        	PairedHit hit = hits.get(i);
        	int leftx = getXPos(hit.leftPos, regionStart, regionEnd, x1, x2);
        	int rightx = getXPos(hit.rightPos, regionStart, regionEnd, x1, x2);
        	int midx = (leftx+rightx)/2;
        	int midy = (int)(((double)(rightx-leftx)/(double)regionWidth) * height*2) + y1;
        	g.setColor(Color.BLACK);
        	g.setStroke(new BasicStroke(1.0f));
        }
	}

}
