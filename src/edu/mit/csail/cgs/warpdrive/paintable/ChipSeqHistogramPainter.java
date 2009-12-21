package edu.mit.csail.cgs.warpdrive.paintable;

import java.io.File;
import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.font.LineMetrics;
import java.util.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.warpdrive.model.ChipSeqHistogramModel;


public class ChipSeqHistogramPainter extends RegionPaintable {

    private ChipSeqHistogramModel model;
    private DynamicAttribute attrib;
    protected static java.util.List configurationFields = null;
    private ChipSeqHistogramProperties props;

    public ChipSeqHistogramPainter(ChipSeqHistogramModel model) {
        super();
        this.model = model;
        props = new ChipSeqHistogramProperties();
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
    }
    public ChipSeqHistogramProperties getProperties() {return props;}
    public void setProperties(ChipSeqHistogramProperties p) {props = p;}
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
    	//long tic =System.currentTimeMillis();
        if (!canPaint()) {
            return;
        }
        if(!model.isReady()) { return; }

        int width = x2 - x1;
        int height = y2 - y1;
        boolean stranded = getProperties().Stranded;
        boolean logscale =false;
        int circlewidth=4;
        
        Map<Integer,Float> plus = model.getPlus(), minus = model.getMinus();

        int maxhits = props.MaxReadCount;
        if (maxhits < 0) {
            double m = 0;
            for (int i : plus.keySet()) {
                if (plus.get(i) > m) {
                    m = plus.get(i);
                }
            }
            for (int i : minus.keySet()) {
                if (minus.get(i) > m) {
                    m = minus.get(i);
                }
            }
            maxhits = (int)m;
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
        if (width / linewidth < model.getRegion().getWidth() / model.getProperties().BinWidth) {
            actualBinWidth = model.getRegion().getWidth() / (width / linewidth);
            combineBins(plus, actualBinWidth);
            combineBins(minus, actualBinWidth);
        }
        int binPixels = width/(model.getRegion().getWidth() / model.getProperties().BinWidth);
        if(binPixels<1)
        	binPixels=1;
        	
        g.setStroke(new BasicStroke((float)linewidth));
        
        // For density probabilities
        if (model.getProperties().GaussianKernelWidth!=0){
        	int totalReadCount = model.getProperties().getTotalReadCount();
            float maxProb = 0;
            for (int i : plus.keySet()) {
                if (plus.get(i) > maxProb) {
                	maxProb = plus.get(i);
                }
            }
            for (int i : minus.keySet()) {
                if (minus.get(i) > maxProb) {
                	maxProb = minus.get(i);
                }
            }
            // re-scale the data
            for (int i : plus.keySet()) {
                plus.put(i, (float)(plus.get(i)/maxProb*totalReadCount*0.5) ); 
            }
            for (int i : minus.keySet()) {
            	minus.put(i, (float)(minus.get(i)/maxProb*totalReadCount*0.5) ); 
            }
            // paint the total read number
            g.setFont(attrib.getLargeLabelFont(width,height));
            g.drawString("Total read: "+totalReadCount,
                         x2-150,
                         y1+15);

        }
        if (stranded) {
            g.setColor(Color.BLUE);
            for (int pos : plus.keySet()) {
                double val = plus.get(pos);
                int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
                int ypix = getYPos(val, 0, maxhits, y1, midpoint, logscale);
                if(model.getProperties().ReadDepthView){
                	System.out.println("+\t"+xpix+"\t"+ypix+"\t"+binPixels+"\t"+actualBinWidth+"\t"+(midpoint-ypix));
                	g.fillRect(xpix, ypix, binPixels, midpoint-ypix);
                }else{
                	g.drawLine(xpix, midpoint, xpix, ypix);
                	g.fillOval(xpix-circlewidth/2, ypix-circlewidth/2, circlewidth, circlewidth);
                }                
            }
            g.setColor(Color.RED);
            for (int pos : minus.keySet()) {
                double val = minus.get(pos);
                int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
                int ypix = midpoint + (y2 - getYPos(val, 0, maxhits, midpoint, y2, logscale));
                if(model.getProperties().ReadDepthView){
                	System.out.println("-\t"+xpix+"\t"+ypix+"\t"+binPixels+"\t"+actualBinWidth+"\t"+(ypix-midpoint));
                	g.fillRect(xpix, midpoint, binPixels, ypix-midpoint);
                }else{
                	g.drawLine(xpix, midpoint, xpix, ypix);
                	g.fillOval(xpix-circlewidth/2, ypix-circlewidth/2, circlewidth, circlewidth);
                }
            }
            g.setColor(Color.black);
            g.drawLine(x1, midpoint, x2, midpoint);
            g.setFont(attrib.getLargeLabelFont(width,height));
            int step = Math.max(1,(int)Math.round(maxhits / 5));
            for (int i = step; i <= Math.ceil(maxhits); i += step) {
                int ypos = getYPos(i, 0, maxhits, y1, midpoint, logscale);
                g.drawString(Integer.toString(i),
                             5,
                             ypos);
                ypos = midpoint + (y2 - getYPos(i, 0, maxhits, midpoint, y2, logscale));
                g.drawString(Integer.toString(i),
                             5,
                             ypos);
            }                    
        } else {
            g.setColor(Color.GRAY);
            for (int pos : plus.keySet()) {
                double val = plus.get(pos);
                if (minus.containsKey(pos)) {
                    val += minus.get(pos);
                }
                int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
                int ypix = getYPos(val, 0, maxhits, y1, y2, logscale);
                if(model.getProperties().ReadDepthView){
                	System.out.println(xpix+"\t"+ypix+"\t"+binPixels+"\t"+actualBinWidth+"\t"+(y2-ypix));
                	g.fillRect(xpix, ypix, binPixels, y2-ypix);
                }else{
                	g.drawLine(xpix, y2, xpix, ypix);
                	g.fillOval(xpix-circlewidth/2, ypix-circlewidth/2, circlewidth, circlewidth);
                }
            }
            for (int pos : minus.keySet()) {
                if (plus.containsKey(pos)) {
                    continue;
                }
                double val = minus.get(pos);
                int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
                int ypix = getYPos(val, 0, maxhits, y1, y2, logscale);
                if(model.getProperties().ReadDepthView){
                	g.fillRect(xpix, ypix, binPixels, y2-ypix);
                }else{
                	g.drawLine(xpix, y2, xpix, ypix);
                	g.fillOval(xpix-circlewidth/2, ypix-circlewidth/2, circlewidth, circlewidth);
                }
            }

            g.setColor(Color.black);
            g.drawLine(x1, y2, x2, y2);
            g.setFont(attrib.getLargeLabelFont(width,height));
            int step = Math.max(1,(int)Math.round(maxhits / 5));
            for (int i = step; i <= Math.ceil(maxhits); i += step) {
                int ypos = getYPos(i, 0, maxhits, y1, y2, logscale);
                g.drawString(Integer.toString(i),
                             5,
                             ypos);
            }                    
        }
        
        g.setStroke(oldStroke);
        if (getProperties().DrawTrackLabel) {
            g.setFont(attrib.getLargeLabelFont(width,height));
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
        }
//        System.err.println((System.currentTimeMillis()-tic)/10*0.01+" sec in painting");
    }

    /**
     * If the binWidth was too small then there aren't enough pixels to show
     * each bin and so they'll be drawn overlapping.  This method
     * resizes the bins into bins of width binWidth
     * by combining skinny bins that map to the same fat bin
     */
    private void combineBins(Map<Integer,Float> map, int binWidth) {
        if (map == null) {
            throw new NullPointerException("null map");
        }
        java.util.List<Integer> keys = new ArrayList<Integer>();
        keys.addAll(map.keySet());
        for (int p : keys) {
            if (p % binWidth == 0) {
                continue;
            } else {
                int newbin = (p / binWidth) * binWidth;
                if (map.containsKey(newbin)) {
                    map.put(newbin, map.get(newbin) + map.get(p));
                } else {
                    map.put(newbin, map.get(p));
                }
                map.remove(p);
            }
        }
    }

}