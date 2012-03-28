package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.io.File;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.warpdrive.model.InteractionLikelihoodModel;

public class InteractionLikelihoodPainter extends RegionPaintable {
	
	private InteractionLikelihoodModel model;
	private InteractionLikelihoodProperties props;
	private DynamicAttribute attrib;
	
	public InteractionLikelihoodPainter(InteractionLikelihoodModel model) {
		super();
		this.model = model;
		props = new InteractionLikelihoodProperties();
		model.addEventListener(this);
		attrib = DynamicAttribute.getGlobalAttributes();
	}
	
	public InteractionLikelihoodProperties getProperties() {
		return props;
	}
	
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
	
	public boolean canPaint() {
		return model.isReady();
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
		int width = x2 - x1;
		int height = y2 - y1;
		Map<Integer,Double> results = model.getResults();
		double minresult = Double.POSITIVE_INFINITY; 
		double maxresult = Double.NEGATIVE_INFINITY;
		for (Integer key : results.keySet()) {
			double tmp = results.get(key);
			if (tmp<minresult) {
				minresult = tmp;
			}
			if (tmp>maxresult) {
				maxresult = tmp;
			}
		}
		System.err.println("minresult: "+minresult+" maxresult: "+maxresult);
		if (maxresult==minresult) {
			minresult = 0;
			maxresult = 1;
		}
		int regionStart = model.getRegion().getStart();
		int regionEnd = model.getRegion().getEnd();
		int midpoint = (y1 + y2) / 2;
		Stroke oldStroke = g.getStroke();
		int linewidth = getProperties().LineWidth;
		if (linewidth < 0) {
			linewidth = width / (model.getRegion().getWidth() / model.getProperties().BinWidth);
		}
		if (linewidth < 1) {
			linewidth = 1;
		}
		int actualBinWidth = model.getProperties().BinWidth;
		int binPixels = width/(model.getRegion().getWidth() / model.getProperties().BinWidth);
		if(binPixels<1)
			binPixels=1;
		g.setStroke(new BasicStroke((float)linewidth));
		
		HashMap<Pair<Integer,Integer>,Double> plotXVals = new HashMap<Pair<Integer,Integer>,Double>();
		g.setColor(Color.GRAY);
		for (int pos : results.keySet()) {
			double val = results.get(pos);
			int xpixleft = getLeftXPos(pos, regionStart, regionEnd, x1, x2);
			int xpixright = getRightXPos(pos, regionStart, regionEnd, x1, x2, actualBinWidth);
			Pair<Integer,Integer> tmppair = new Pair<Integer,Integer>(xpixleft,xpixright);
			if(!plotXVals.containsKey(tmppair) || plotXVals.get(tmppair)<val){
				plotXVals.put(tmppair,val);
			}
		}
		for(Pair<Integer,Integer> tmppair : plotXVals.keySet()){
			double val = plotXVals.get(tmppair); 
			int ypix = getYPos(val, minresult, maxresult, y1, y2, false);
			g.fillRect(tmppair.car(), ypix, tmppair.cdr()-tmppair.car(), y2-ypix);
		}
		System.err.println();
		//Line & trimmings
		g.setColor(Color.black);
		g.drawLine(x1, y2, x2, y2);
		g.setFont(attrib.getLargeLabelFont(width,height));
		double step = (maxresult-minresult)/5d;
		for (double i = step; i <= maxresult; i += step) {
			int ypos = getYPos(i, minresult, maxresult, y1, y2, false);
			g.drawString(Double.toString(i),5,ypos);
		}
		
		g.setColor(Color.red);
		Region anchor = Region.fromString(model.getRegion().getGenome(), model.getProperties().Anchor);
		int anchorx1 = getXPos(anchor.getStart(), regionStart, regionEnd, x1, x2);
		g.drawLine(anchorx1, y1, anchorx1, y2);
		anchorx1 = getXPos(anchor.getEnd(), regionStart, regionEnd, x1, x2);
		g.drawLine(anchorx1, y1, anchorx1, y2);
		
		g.setStroke(oldStroke);
		if (getProperties().DrawTrackLabel) {
			g.setFont(attrib.getLargeLabelFont(width,height));
			g.setColor(Color.BLACK);
			g.drawString(getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}if(getProperties().DrawBinSize) {
			g.setFont(attrib.getPointLabelFont(width,height));
			g.setColor(Color.GRAY);
			String binString = new String("(bin size: "+actualBinWidth+"bp)");
			g.drawString(binString,x2 - g.getFontMetrics().stringWidth(binString) - g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}
	}
	
	public int getLeftXPos(int pos, int start, int end, int leftx, int rightx) {
        if (pos < start) {return leftx;}
        if (pos > end) {return rightx;}
        return (int)((((float)(pos - start))/((float)(end - start))) * (rightx - leftx) + leftx);
    }
	
	public int getRightXPos(int pos, int start, int end, int leftx, int rightx, int binwidth) {
		if (pos > end) {return rightx;}
		if (pos < start) {
        	while (pos<start) {
        		pos += binwidth;
        	}
        } else {
        	pos += binwidth;
        }
        return (int)((((float)(pos - start))/((float)(end - start))) * (rightx - leftx) + leftx);
    }

}
