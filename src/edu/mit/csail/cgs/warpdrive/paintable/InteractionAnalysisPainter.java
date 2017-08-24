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

	static final private Color[] arcColors = {new Color(0, 0, 255, 127), new Color(0, 255, 255, 127), 
			new Color(0, 255, 128, 127), new Color(255, 128, 0, 127), new Color(255, 0, 0, 127)};
	static final private Color[] textColors = {new Color(0, 0, 255, 255), new Color(0, 255, 255, 255), 
			new Color(0, 255, 128, 255), new Color(255, 128, 0, 255), new Color(255, 0, 0, 255)};
	static final private Color lightGrey = new Color(0.4f, 0.4f, 0.4f, 0.2f);
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
		String chrom = model.getRegion().getChrom();
		int regionWidth = model.getRegion().getWidth();
		int linewidth = Math.max(getProperties().LineWidth,1);
		Stroke oldStroke = g.getStroke();
		g.setStroke(new BasicStroke((float)linewidth));

		Map<Pair<Point,Point>,Float> interactions = model.getInteractions();
		Map<Pair<Point,Point>, Pair<Region,Region>> interactionAnchors = model.getInteractionAnchors();
		
		if (getProperties().DrawTrackLabel) {
			g.setFont(attrib.getLargeLabelFont(width,height));
			g.setColor(Color.BLACK);
			g.drawString("Analysis " +getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}

		g.setStroke(new BasicStroke(1.0f));
		TreeMap<Point, Integer> leftOutRanges = new TreeMap<Point, Integer>();
		TreeMap<Point, Integer> rightOutRanges = new TreeMap<Point, Integer>();
		int cutoff = getProperties().ReadPairCountCutoff;
		boolean toShowCount = getProperties().DisplayReadPairCount;
		boolean toShowOutOfRangeIntereactions = getProperties().ShowOutOfRangeIntereactions;
		boolean toShowAnchors = getProperties().ShowAnchors;
		for (Pair<Point,Point> pair : interactions.keySet()) {
			Point leftPoint = pair.car();
			Point rightPoint = pair.cdr();
			if (leftPoint.equals(rightPoint))
				continue;
			if (!leftPoint.getChrom().equals(chrom))
				continue;
			float count = interactions.get(pair);
			if (count<cutoff && count!=1)
				continue;
			float curvewidth = Math.min(30, (float)Math.sqrt((double) count));
			g.setStroke(new BasicStroke(curvewidth, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));
			int leftCoord = leftPoint.getLocation();
			int rightCoord = rightPoint.getLocation();
			if (leftCoord<regionStart){
				// if out of range, record the position that is in range
				if (!leftOutRanges.containsKey(rightPoint))
					leftOutRanges.put(rightPoint, 1);
				else
					leftOutRanges.put(rightPoint, leftOutRanges.get(rightPoint)+1);
				if (!toShowOutOfRangeIntereactions)
					continue; // skip
			}
			if (rightCoord>regionEnd){
				// if out of range, record the position that is in range
				if (!rightOutRanges.containsKey(leftPoint))
					rightOutRanges.put(leftPoint, 1);
				else
					rightOutRanges.put(leftPoint, rightOutRanges.get(leftPoint)+1);
				if (!toShowOutOfRangeIntereactions)
					continue; // skip
			}
			int leftx = getXPosExt(leftCoord, regionStart, regionEnd, x1, x2);
			int rightx = getXPosExt(rightCoord, regionStart, regionEnd, x1, x2);
			int midx = (leftx+rightx)/2;
			int midy = y2 - (int)(((double)(rightx-leftx)/(double)width) * 2*height);
			int colorIdx = 0;
			if (count==1){
				g.setColor(lightGrey);
				g.setStroke(new BasicStroke(1/15f));
			}
			else{
				colorIdx = Math.round(count/10);
				if (colorIdx > arcColors.length-1)
					colorIdx = arcColors.length-1;
				g.setColor(arcColors[colorIdx]);
				// draw anchors
				if (toShowAnchors){
					Pair<Region, Region> anchors = interactionAnchors.get(pair);
					Region leftAnchor = anchors.car();
					Region rightAnchor = anchors.cdr();
					int leftS = getXPosExt(leftAnchor.getStart(), regionStart, regionEnd, x1, x2);
					int leftE = getXPosExt(leftAnchor.getEnd(), regionStart, regionEnd, x1, x2);
					int rightS = getXPosExt(rightAnchor.getStart(), regionStart, regionEnd, x1, x2);
					int rightE = getXPosExt(rightAnchor.getEnd(), regionStart, regionEnd, x1, x2);
					g.fillRect(leftS, y2, leftE-leftS, 3);
					g.fillRect(rightS, y2, rightE-rightS, 3);
				}
			}
			QuadCurve2D loop = new QuadCurve2D.Float(leftx, y2, midx, midy, rightx, y2);
			g.draw(loop);
			if (toShowCount && count!=1){
				g.setColor(textColors[colorIdx]);
				String countStr = count==Math.round(count) ? Math.round(count)+"" : count+"";
				g.drawString(countStr, midx, y2 - (int)(((double)(rightx-leftx)/(double)width) * height));
			}
		}
		// list out-of-range interactions
		if (toShowOutOfRangeIntereactions) {
			g.setFont(attrib.getPointLabelFont(width,height));
			int fontSize = g.getFont().getSize();
			g.setColor(Color.RED);
			g.drawString("<", x1, y1 + (int)(fontSize*2.5));
			for (Point p: leftOutRanges.keySet()){
				g.drawString(leftOutRanges.get(p)+".", 
						getXPosExt(p.getLocation(), regionStart, regionEnd, x1, x2), 
						y1 + (int)(fontSize*2.5));
			}
			g.setColor(Color.BLACK);
			g.drawString(">", x2-fontSize, y1 + fontSize*4);
			for (Point p: rightOutRanges.keySet()){
				g.drawString(rightOutRanges.get(p)+".", 
						getXPosExt(p.getLocation(), regionStart, regionEnd, x1, x2),
						y1 + fontSize*4);
			}		
		}
		g.setStroke(oldStroke);
	}
	/** override to return real scaled out_of_range positions */
    private int getXPosExt(int pos, int start, int end, int leftx, int rightx) {
        return (int)((((float)(pos - start))/((float)(end - start))) * (rightx - leftx) + leftx);
    }


	public InteractionAnalysisProperties getProperties() {
		return props;
	}

}
