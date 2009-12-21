package edu.mit.csail.cgs.warpdrive.paintable;

import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.awt.font.LineMetrics;
import java.io.File;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.GenericExperiment;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.warpdrive.model.ChipChipDataModel;
import edu.mit.csail.cgs.warpdrive.model.ChipChipScaleModel;
import edu.mit.csail.cgs.warpdrive.model.Model;
import edu.mit.csail.cgs.warpdrive.model.ModelProperties;
import edu.mit.csail.cgs.ewok.nouns.*;

/* paints ChIP-Chip data across a genomic region. */

public class TimChipChipPainter extends RegionPaintable  {

    public static int CHIP = 1, EXPRESSION = 2, CGH = 3, RULER = 4;
    public static ColorSet colors = new ColorSet();
    private GenericExperiment data;
    private ChipChipDataModel model;
    private ChipChipScaleModel scale;
    private DynamicAttribute attrib;
    private int circleRadius, circleWidth, prevmaxcount;
    private static java.util.List configurationFields = null;
    private ExtendedChipChipProperties props;
    private Color color; // this is basically a temporary variable shared by several methods

    public TimChipChipPainter (GenericExperiment data, ChipChipDataModel model) {
        super();
        this.data = data;
        this.model = model;
        props = new ExtendedChipChipProperties();
        props.Color = colors.getColor(this.toString());
        model.addEventListener(this);
        attrib = DynamicAttribute.getGlobalAttributes();
        prevmaxcount = -1;
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }
    
    public void setScaleModel(ChipChipScaleModel s) {
        scale = s;
    }
    public ExtendedChipChipProperties getProperties() {return props;}
    public void savePropsInDir(File dir) {
        super.savePropsInDir(dir);
        ModelProperties mp = model.getProperties();
        File saveto = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + mp.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        mp.saveToFile(saveto);
    }
    public void loadPropsInDir(File dir) {
        super.loadPropsInDir(dir);
        ModelProperties mp = model.getProperties();
        File loadfrom = new File(dir + System.getProperty("file.separator") + (getLabel() + "." + mp.defaultName()).replaceAll(System.getProperty("file.separator"),"_"));
        mp.saveToFile(loadfrom);
    }
    public ChipChipScaleModel getScaleModel() {return scale;}

    /* style refers to one of the RULER, EXPRESSION, CGH, etc constants
     */
    public void setChannelStyle(int style) {
        synchronized(props) {
            if (style == RULER || style == CGH || style == EXPRESSION) {
                getProperties().DrawCy5 = true;
                getProperties().DrawCy3 = true;
                getProperties().DrawStrandedRatio = true;
            }
            if (style == CGH) {
                getProperties().MaxRatio = 4.0;
                getProperties().MaxIntensity = 20000.0;
            }
            if (style == EXPRESSION) {
                getProperties().MaxIntensity = 15000.0;
            }
        }
    }

    public java.util.List<String> configurationKeyOrder() {return configurationFields;}

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        int w = x2 - x1, h = y2 - y1;
        Region region = model.getRegion();
        int regionWidth = region.getEnd() - region.getStart() + 1;
        int rs = region.getStart(), re = region.getEnd();
        int smoothpoints = getProperties().AverageRatiosAcrossNPoints;
        if (smoothpoints < 1) {smoothpoints = 1;}

    	//--------------------------------------------------

        boolean logscale = getProperties().RatiosOnLogScale;
        boolean logscaleintensities = getProperties().IntensitiesOnLogScale;
        int maxLineDistance = getProperties().MaximumLineDistance;
        maxLineDistance = (int)(maxLineDistance * ((float)w) / ((float)regionWidth));
        Color lineColor = getProperties().Color;
        ArrayList<Color> repColors = new ArrayList<Color>();
        double maxValue;
        if (getProperties().MaxRatio > 0) {
            maxValue = getProperties().MaxRatio;
        } else {
            maxValue = scale.getMaxVal();
        }
        double minValue = 0.0;        
        if (getProperties().RatiosOnLogScale) {
            minValue = 1.0 / maxValue;
        }
                
        if(!model.isReady()) { 
            g.setColor(Color.black);
            g.drawString("Loading... Please Wait.", x1 + w / 2, y1 + h / 2);
            return;
        }
        
        Graphics2D g2 = (Graphics2D)g;
        Stroke oldStroke = g2.getStroke();

        int ppx = -1, ppy = -1;
        g2.setColor(lineColor);
        
        GenericExperiment expt = model.getGenericExperiment();

        double drawingscale = 1.0 / Math.log(Math.max((double)expt.getCount() * 5, 1.0));

        int lineWidth = attrib.getLineWidth(w, h, drawingscale);
        int propLineWidth = getProperties().MinLineWidth;
        lineWidth = Math.max(lineWidth, propLineWidth);
        
        circleRadius = attrib.getPointWidth(w, h, drawingscale);
        int minCircleRadius = getProperties().MinCircleRadius;
        circleRadius = Math.max(minCircleRadius, circleRadius);
        int circleLineWidth = circleRadius;        
        circleWidth = circleRadius * 2;
        
        Stroke lineStroke = new BasicStroke((float)lineWidth);
        Stroke circleStroke = new BasicStroke((float)circleLineWidth);        
        g2.setStroke(lineStroke);
        
        Font oldFont = g.getFont();
        Font newFont = attrib.getRegionLabelFont(w, h);
        int fontSize = getProperties().FontSize;
        if(fontSize > 0 && fontSize != newFont.getSize()) { 
        	newFont = new Font(newFont.getName(), newFont.getStyle(), fontSize);
        }
        g.setFont(newFont);
        
        if (getProperties().DrawRatios) {
            boolean average = getProperties().AverageAcrossReplicates;
            boolean bystrand = getProperties().DrawStrandedRatio;

            int[][][] locs = new int[expt.getCount()][][];
            char[][] strand = null;
            if (bystrand && expt instanceof ChipChipData) {
                ChipChipData ccd = (ChipChipData)expt;
                strand = new char[ccd.getCount()][];
                for (int i = 0; i < ccd.getCount(); i++) {
                    strand[i] = new char[ccd.getReplicates(i)];
                    for (int j = 0; j < ccd.getReplicates(i); j++) {
                        strand[i][j] = ccd.getStrand(i,j);
                    }

                }

            }

            if (repColors.size() == 0) {
                repColors.add(lineColor.brighter());
            }
            if (average) {
                lineStroke = new BasicStroke((float)lineWidth * 2);
            }
            g.setStroke(lineStroke);
            for(int i = 0; i < expt.getCount(); i++) {
                int px = getXPos(expt.getPos(i),
                                 rs,re,x1,x2);
                int reps = expt.getReplicates(i);
                char thisstrand = bystrand ? strand[i][0] : ' ';
                int lowerbound = (int)Math.max(0, i - Math.floor(smoothpoints / 2.0));
                int upperbound = (int)Math.min(locs.length - 1, i + Math.ceil(smoothpoints / 2.0) - 1);
                if (average) { // average across replicates
                    locs[i] = new int[1][2];                   
                    double sum = 0;
                    double count = 0;
                    for (int k = lowerbound; k <= upperbound; k++) {
                        reps = expt.getReplicates(k);
                        for(int j = 0; j < reps; j++) { 
                            if (bystrand && strand[k][j] != thisstrand) {continue;}
                            double v = expt.getValue(k,j);
                            if (!Double.isNaN(v)) {
                                sum += v;
                                count++;
                            }
                        }
                    }
                    double value = sum / count;
                    int py = getYPos(value,
                                     minValue,maxValue,
                                     y1,y2, logscale);                        
                    if (Double.isNaN(value)) {
                        py = -1;
                    }
                    locs[i][0][0] = px; locs[i][0][1] = py;                        
                } else {
                    locs[i] = new int[expt.getReplicates(i)][2];
                    for(int j = 0; j < reps; j++) { 
                        double sum = 0;
                        double count = 0;
                        thisstrand = bystrand ? strand[i][j] : ' ';
                        for (int k = lowerbound; k <= upperbound; k++) {
                            if (j >= expt.getReplicates(k)) {
                                continue;
                            }
                            if (bystrand && strand[k][j] != thisstrand) {continue;}
                            double val = expt.getValue(k,j);
                            if (!Double.isNaN(val)) {
                                sum  += val;
                                count++;
                            }
                        }
                        double value = sum / count;
                        int py = getYPos(value,
                                         minValue,maxValue,
                                         y1,y2, logscale);
                        if (Double.isNaN(value)) {
                            py = -1;
                        }
                        
                        locs[i][j][0] = px; locs[i][j][1] = py;
                        if (j >= repColors.size()) {
                            repColors.add(repColors.get(j-1).darker());
                        }
                    }                        
                }
            }
            if (bystrand && locs.length != strand.length) {
                throw new RuntimeException(String.format("locs length is %d but strand length is %d",
                                                         locs.length, strand.length));
            }
            for (int i = 0; i < locs.length - 1; i++) {
                for (int j = 0; j < locs[i].length; j++) {
                    if (locs[i][j][1] < 0) {
                        continue;
                    }
                    g.setColor(repColors.get(j));                    
                    int k = i;
                    while (k < i + 10 && k < locs.length - 1) {
                        k++;
                        if (locs[i].length != locs[k].length) {
                            continue;
                        }
                        if (bystrand) {
                            if (strand[k].length != strand[i].length) {
                                continue;
                            }
                            if (strand[k][j] != strand[i][j]) {
                                continue;
                            }
                        }
                        if ((Math.abs(locs[i][j][0] - locs[k][j][0]) < maxLineDistance) &&                            
                            (locs[k][j][1] >= 0)) {
                            g2.drawLine(locs[k][j][0], locs[k][j][1], locs[i][j][0], locs[i][j][1]);
                            break;
                        } 
                    }
                }
            }

            g2.setStroke(circleStroke);
            for(int i = 0; i < locs.length; i++) {
                for(int j = 0; j < locs[i].length; j++) { 
                    color = repColors.get(j);
                    if (locs[i][j][1] < 0) {
                        g.fillOval(locs[i][j][0] - circleRadius,
                                   y2 - circleRadius,
                                   circleWidth,
                                   circleWidth);
                    } else {
                        paintDataPointAt(g,locs[i][j][0], locs[i][j][1],i,j);
                    }
                }
            }
        }
        boolean drawip = getProperties().DrawCy5, drawwce = getProperties().DrawCy3;
        if ((model.getGenericExperiment() instanceof ChipChipData) && (drawip || drawwce)) {

            int smoothintensitypoints = getProperties().AverageIntensitiesAcrossNPoints;
            ChipChipData data = (ChipChipData)model.getGenericExperiment();
            double maxval = getProperties().MaxIntensity;
            double minval = getProperties().MinIntensity;
            if (maxval < 0 || maxval < minval) {
                for (int i = 0; i < data.getCount(); i++) {
                    for (int j = 0; j < data.getReplicates(i); j++) {
                        if (data.getWCE(i,j) > maxval && drawwce) {
                            maxval = data.getWCE(i,j);
                        }
                        if (data.getIP(i,j) > maxval && drawip) {
                            maxval = data.getIP(i,j);
                        }
                    }
                }
            }

            g.setStroke(lineStroke);
            /* draw the scale for intensities */
            g.setFont(DynamicAttribute.getGlobalAttributes().getRegionLabelFont(w,h));
            FontMetrics fontmetrics = g.getFontMetrics();
            int marks = h / (4 * fontmetrics.getHeight());
            double mult = Math.max(Math.exp(Math.log(maxval / minval) / marks), 2);
            g.setColor(Color.BLACK);
            for (double d = minval; d <= maxval; d *= mult) {
                int y = getYPos(d,minval,maxval,y1,y2, logscaleintensities);
                String l = Integer.toString((int)d);
                int lwidth = fontmetrics.charsWidth(l.toCharArray(),0,l.length());                
                g.drawLine(x2-10,y,x2,y);
                g.drawString(l,x2-10 - lwidth,y + fontmetrics.getHeight()/2);
            }
            
            for (int i = 0; i < data.getCount(); i++) {
                int px = getXPos(data.getPos(i),
                                 rs,re,x1,x2);
                for (int j = 0; j < data.getReplicates(i); j++) {
                    char strand = data.getStrand(i,j);
                    int lowerbound = (int)Math.max(0, i - Math.floor(smoothintensitypoints / 2.0));
                    int upperbound = (int)Math.min(data.getCount() - 1, i + Math.ceil(smoothintensitypoints / 2.0) - 1);
                    double ipval = 0, wceval = 0;
                    int count = 0;
                    for (int k = lowerbound; k <= upperbound; k++) {
                        if (j >= expt.getReplicates(k)) {continue;}
                        if (!Double.isNaN(data.getIP(k,j)) && !Double.isNaN(data.getWCE(k,j))) {
                            ipval += data.getIP(k,j);
                            wceval += data.getWCE(k,j);
                            count++;
                        }
                    }
                    ipval /= count;
                    wceval /= count;

                    if (drawip && !Double.isInfinite(ipval) && !Double.isNaN(ipval)) {
                        g.setColor(Color.RED);
                        if (Double.isNaN(ipval)) {
                            ipval = minval;
                            g.setColor(Color.GRAY);
                            strand = 'n';
                        }
                        int ippy = getYPos(ipval,
                                           minval,maxval,
                                           y1,y2,
                                           logscaleintensities);  

                        if (strand == '+') {
                            plusAt(px-circleRadius,ippy-circleRadius,circleWidth,circleWidth,g);
                        } else if (strand == '-') {
                            minusAt(px-circleRadius,ippy-circleRadius,circleWidth,circleWidth,g);
                        } else {
                            ovalAt(px-circleRadius,ippy-circleRadius,circleWidth,circleWidth,g);
                        }
                    }
                    if (drawwce && !Double.isInfinite(wceval) && !Double.isNaN(wceval)) {
                        g.setColor(Color.GREEN);
                        if (Double.isNaN(wceval)) {
                            wceval = minval;
                            g.setColor(Color.GRAY);
                            strand = 'n';
                        }
                        int wcepy = getYPos(wceval,
                                            minval,maxval,
                                            y1,y2, logscaleintensities);
                        
                        if (strand == '+') {
                            plusAt(px-circleRadius,wcepy-circleRadius,circleWidth,circleWidth,g);
                        } else if (strand == '-') {
                            minusAt(px-circleRadius,wcepy-circleRadius,circleWidth,circleWidth,g);
                        } else {
                            ovalAt(px-circleRadius,wcepy-circleRadius,circleWidth,circleWidth,g);
                        }
                    }
                }
            }
        }
        g2.setStroke(oldStroke);        
        if (getProperties().DrawTrackLabel) {        	
            g.setFont(attrib.getLargeLabelFont(w,h));
            g.setColor(Color.BLACK);
            boolean onRight = getProperties().DrawLabelOnRight;
            
            if(onRight) { 
            	FontRenderContext frc = g2.getFontRenderContext();
            	TextLayout layout = new TextLayout(getLabel(), g.getFont(), frc);
            	int width = (int)Math.ceil(layout.getBounds().getWidth());
            	g.drawString(getLabel(),x2 - width - 2, y1 + g.getFont().getSize() * 2);
            } else {                    
            	g.drawString(getLabel(),x1,y1 + g.getFont().getSize() * 2);
            }
        }

        g.setFont(oldFont);
    }

    private void plusAt(int x, int y, int w, int h, Graphics2D g) {
        g.drawLine(x-w,y,x+w,y);
        g.drawLine(x,y-h,x,y+h);
    }
    private void minusAt(int x, int y, int w, int h, Graphics2D g) {
        g.drawLine(x-w,y,x+w,y);
        g.drawLine(x-w,y-1,x+w,y-1);
        g.drawLine(x-w,y+1,x+w,y+1);
    }
    private void ovalAt(int x, int y, int w, int h, Graphics2D g) {
        g.fillOval(x,y,w,h);
    }

    public void paintDataPointAt(Graphics2D g, int x,int y,int i,int j) {        
        g.setColor(color);
        g.fillOval(x - circleRadius, y - circleRadius, circleWidth, circleWidth);
        g.setColor(Color.black);
        g.drawOval(x - circleRadius, y - circleRadius, circleWidth, circleWidth);
    }

    public synchronized void eventRegistered(EventObject e) {        
        if ((e.getSource() == model) || (e.getSource() == scale) &&
            model.isReady() &&
            ((scale == null) || scale.isReady())) {
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

