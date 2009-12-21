package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.font.LineMetrics;
import java.util.*;
import java.text.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.*;

public class ExpressionMeasurementPainter extends RegionPaintable  {
	
	private static NumberFormat nf;
    private PaintableProperties props;
	
	static {
		nf = DecimalFormat.getInstance();
		nf.setMaximumFractionDigits(2);
		nf.setMinimumFractionDigits(2);
	}

    private ExpressionProbeModel model;
    private DynamicAttribute attrib;
    private static java.util.List configurationFields = null;

    public ExpressionMeasurementPainter (ExpressionProbeModel m) {
        super();
        model = m;
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
        props = new PaintableProperties();
    }
    
    public PaintableProperties getProperties() {
        return props;
    }

    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public java.util.List<String> configurationKeyOrder() {return configurationFields;}

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
    	int w = x2 - x1, h = y2 - y1;
        
        if(!model.isReady()) { 
            g.setColor(Color.black);
            g.drawString("Loading... Please Wait.", x1 + w / 2, y1 + h / 2);
            return;
        }
        
        Iterator<LocatedExprMeasurement> itr = model.getResults();
        LinkedList<LocatedExprMeasurement> ems = new LinkedList<LocatedExprMeasurement>();
        double maxValue = 0.0;
        while(itr.hasNext()) {
        	LocatedExprMeasurement em = itr.next();
        	ems.addLast(em);
        	maxValue = Math.max(maxValue, em.getValue());
        }
        
        System.out.println("Expression: " + ems.size() + " values.");
        
        Rectangle rect = new Rectangle(x1, y1, w, h);
        Region r = model.getRegion();
        
        for(LocatedExprMeasurement em : ems) { 
        	paintMeasurement(g, r, rect, em, maxValue);
        }
    }
    
    private void paintMeasurement(Graphics2D g, Region r, Rectangle bounds, 
    		LocatedExprMeasurement em, double maxValue) { 
    	
    	Region emr = em.getRegion();

    	if(emr.getEnd() < r.getStart() || emr.getStart() > r.getEnd()) { 
    		System.err.println("Error: " + em);
    		return;
    	}
    	
    	double val = em.getValue();
    	double valf = val / maxValue;
    	int barh = Math.max(1, (int)Math.round(valf * (double)bounds.height));
    	
    	int rw = r.getEnd() - r.getStart() + 1;
    	
    	double sf = (double)(emr.getStart()-r.getStart()) / (double)rw;
    	double ef = (double)(emr.getEnd() - r.getStart()) / (double)rw;
    	
    	int barstart = bounds.x + (int)Math.round(sf * (double)bounds.width);
    	int barend = bounds.x + (int)Math.round(ef * (double)bounds.width);
    	
    	int barwidth = Math.max(1, barend - barstart + 1);
    	
    	int bary = bounds.y + bounds.height - barh;
		int radius = 2;
		int diam = 2 * radius;

        int barx = (barstart+barend)/2;
        g.setColor(Color.pink);
        g.drawLine(barx, bounds.y + bounds.height, barx, bary);
        g.setColor(Color.white);
        g.fillOval(barx-radius, bary-radius, diam, diam);
        g.setColor(Color.black);
        g.drawOval(barx-radius, bary-radius, diam, diam);
        
        /*
    	if(barwidth > 1) { 
    		g.setColor(Color.pink);
    		g.fillRect(barstart, bary, barwidth, barh);
    		g.setColor(Color.black);
    		g.drawRect(barstart, bary, barwidth, barh);
    	} else { 
    		g.setColor(Color.black);
    		g.drawLine(barstart, bounds.y + bounds.height, barstart, bary);
    		g.drawOval(barstart-radius, bary-radius, diam, diam);
    	}
        */
    	
    	//g.drawString(nf.format(em.getValue()), barstart, bary - diam);
    }

    public synchronized void eventRegistered(EventObject e) {        
        if ((e.getSource() == model) &&
            model.isReady()) {
        	
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
}

