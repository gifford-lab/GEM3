package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.Graphics2D;
import java.awt.Stroke;
import java.io.File;
import java.util.*;

import edu.mit.csail.cgs.utils.Listener;
import edu.mit.csail.cgs.warpdrive.model.InteractionLikelihoodModel;

public class InteractionLikelihoodPainter extends RegionPaintable {
	
	private InteractionLikelihoodModel model;
	private InteractionLikelihoodProperties props;
	
	public InteractionLikelihoodPainter(InteractionLikelihoodModel model) {
		super();
		this.model = model;
		props = new InteractionLikelihoodProperties();
		model.addEventListener(this);
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
	}

}
