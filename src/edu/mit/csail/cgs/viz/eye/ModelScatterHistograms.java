/*
 * Author: tdanford
 * Date: Nov 14, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.NoninvertibleTransformException;
import java.util.Iterator;

import edu.mit.csail.cgs.utils.models.Model;

public class ModelScatterHistograms extends ContainerModelPaintable {
	
	private ModelScatter scatter;
	private ModelHistogram xHist, yHist;
	
	public ModelScatterHistograms(String xName, String yName) {
		scatter = new ModelScatter(xName, yName);
		xHist = new ModelHistogram(xName);
		yHist = new ModelHistogram(yName);
		
		addModelPaintable(scatter);
		addModelPaintable(xHist);
		addModelPaintable(yHist);
		
		xHist.setProperty(ModelHistogram.colorKey, Color.blue);
		yHist.setProperty(ModelHistogram.colorKey, Color.blue);

		scatter.setProperty(ModelScatter.xScaleKey, scatter.getProperty(ModelScatter.yScaleKey));
		scatter.synchronizeProperty(ModelScatter.xScaleKey, xHist, ModelHistogram.xScaleKey);
		scatter.synchronizeProperty(ModelScatter.yScaleKey, yHist, ModelHistogram.xScaleKey);
		
		xHist.setProperty(ModelHistogram.binsKey, 50);
		yHist.setProperty(ModelHistogram.binsKey, 50);
	}
	
	public ModelScatter getModelScatter() { return scatter; }
	public ModelHistogram getXHistogram() { return xHist; }
	public ModelHistogram getYHistogram() { return yHist; }

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		int w4 = w/5, h4 = h/5;
		
		scatter.paintItem(g, x1+w4, y1, x2, y2-h4);
		xHist.paintItem(g, x1+w4, y2-h4, x2, y2);
		
		Graphics2D g2 = (Graphics2D)g;
		
		//int dx = 0, dy = y2-h4;
		//double theta = -Math.PI/2.0;

		int dx = x1+w4, dy = y2-h4;
		double theta = -3.0 * Math.PI/2.0;
		AffineTransform transform = null, inverse = null;
		try {
			transform = new AffineTransform(-1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
			inverse = transform.createInverse();
		} catch (NoninvertibleTransformException e) {
			e.printStackTrace();
			transform = inverse = null;
		}
		
		g2.translate(dx, dy);
		g2.rotate(theta);
		if(transform != null) { g2.transform(transform); } 
		
		yHist.paintItem(g, 0, 0, h-h4, w4);
		
		if(inverse != null) { g2.transform(inverse); }
		g2.rotate(-theta);
		g2.translate(-dx, -dy);
	}
	
	public void rebin() { 
		xHist.rebin();
		yHist.rebin();
	}
	
	public void addModel(Model m) { 
		super.addModel(m);
	}
	
	public void addModels(Iterator<? extends Model> ms) { 
		super.addModels(ms);
		rebin();
	}
}
