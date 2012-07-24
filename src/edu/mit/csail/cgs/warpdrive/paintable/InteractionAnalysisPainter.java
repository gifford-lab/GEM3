package edu.mit.csail.cgs.warpdrive.paintable;

import edu.mit.csail.cgs.datasets.general.*;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.QuadCurve2D;
import java.util.*;

import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.warpdrive.model.InteractionAnalysisModel;

public class InteractionAnalysisPainter extends RegionPaintable {

	private InteractionAnalysisModel model;
	private InteractionAnalysisProperties props;
	private DynamicAttribute attrib;

	public InteractionAnalysisPainter(InteractionAnalysisModel model) {
		super();
		this.model = model;
		props = new InteractionAnalysisProperties();
		model.addEventListener(this);
		attrib = DynamicAttribute.getGlobalAttributes();
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
		int halfy = y2 - (height/2);
		int regionStart = model.getRegion().getStart();
		int regionEnd = model.getRegion().getEnd();
		int regionWidth = model.getRegion().getWidth();
		int linewidth = Math.max(getProperties().LineWidth,1);
		Stroke oldStroke = g.getStroke();
		g.setStroke(new BasicStroke((float)linewidth));

		Map<Point,Float> events = model.getEvents();
		Map<Pair<Point,Point>,Float> interactions = model.getInteractions();

		if (getProperties().DrawTrackLabel) {
			g.setFont(attrib.getLargeLabelFont(width,height));
			g.setColor(Color.BLACK);
			g.drawString("Analysis " +getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}

		int h = height;
		float maxweight = 0;
		g.setStroke(new BasicStroke(1.0f));
		for (Pair<Point,Point> pair : interactions.keySet()) {
			if (!pair.car().equals(pair.cdr())) {
				float count = interactions.get(pair);
				float curvewidth = Math.min(30, count);
				g.setStroke(new BasicStroke(curvewidth, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));
				int leftx = getXPos(pair.car().getLocation(), regionStart, regionEnd, x1, x2);
				int rightx = getXPos(pair.cdr().getLocation(), regionStart, regionEnd, x1, x2);
				int midx = (leftx+rightx)/2;
				int midy = y2 - (int)(((double)(rightx-leftx)/(double)width) * 2*height);
				g.setColor(new Color(0, 0, 255, 127));
				QuadCurve2D loop = new QuadCurve2D.Float(leftx, y2, midx, midy, rightx, y2);
				g.draw(loop);
				g.setColor(Color.black);
				g.drawString(count+"", midx, y2 - (int)(((double)(rightx-leftx)/(double)width) * height));
			}
		}
		g.setStroke(oldStroke);
	}

	public InteractionAnalysisProperties getProperties() {
		return props;
	}

}
