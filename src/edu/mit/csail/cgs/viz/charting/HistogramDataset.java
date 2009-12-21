/*
 * Created on Aug 16, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.viz.charting;

import java.util.HashMap;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.IntervalXYDataset;

public class HistogramDataset implements IntervalXYDataset {
    
    private static Logger logger = Logger.getLogger("edu.mit.csail.cgs.viz.charting.HistogramDataset");
    
    private double min, max;
    private double binWidth;
    private int bins;
    
    private Vector<String> seriesNames;
    private double[] xstarts, xends;
    private Map<String,int[]> counts;
    private Map<String,Integer> seriesCounts;
    
    public HistogramDataset(double min, double max, int bins) { 
        this.min = min; this.max = max;
        this.bins = bins;
        binWidth = (max-min)/(double)bins;
        
        seriesNames = new Vector<String>();
        xstarts = new double[bins];
        xends = new double[bins];
        counts = new HashMap<String,int[]>();
        seriesCounts = new HashMap<String,Integer>();
        
        double bstart = min, bend = bstart + binWidth;
        for(int i = 0; i < bins; i++, bstart += binWidth, bend += binWidth) { 
            xstarts[i] = bstart;
            xends[i] = bend;
        }
    }
    
    public int addSeries(String name) {
        if(seriesNames.contains(name)) { throw new IllegalArgumentException(name); }
        
        seriesNames.add(name);
        seriesCounts.put(name, 0);
        int[] array = new int[bins];
        for(int i = 0; i < bins; i++) { array[i] = 0; }
        counts.put(name, array);
        
        return seriesNames.size()-1;
    }
    
    public void addValue(int series, double value) { 
        int bin = findBin(value);
        String seriesName = seriesNames.get(series);
        counts.get(seriesName)[bin] += 1;
        seriesCounts.put(seriesName, seriesCounts.get(seriesName) + 1);
    }
    
    private int findBin(double val) { 
        if(val < min || val > max) { 
            return -1;
        }
        int bin = (int)Math.floor((val-min) / binWidth);
        return Math.min(bin, bins-1);
    }

    public Number getEndX(int s, int i) {
        return xends[i];
    }

    public double getEndXValue(int s, int i) {
        return xends[i];
    }

    public Number getEndY(int s, int i) {
        return counts.get(seriesNames.get(s))[i];
    }

    public double getEndYValue(int s, int i) {
        return counts.get(seriesNames.get(s))[i];
    }

    public Number getStartX(int s, int i) {
        return xstarts[i];
    }

    public double getStartXValue(int s, int i) {
        return xstarts[i];
    }

    public Number getStartY(int s, int i) {
        return 0.0;
    }

    public double getStartYValue(int arg0, int arg1) {
        return 0.0;
    }

    public DomainOrder getDomainOrder() {
        return DomainOrder.NONE;
    }

    public int getItemCount(int s) {
        return bins;
    }

    public Number getX(int s, int i) {
        return (xstarts[i] + xends[i])/2.0; 
    }

    public double getXValue(int s, int i) {
        return (xstarts[i] + xends[i]) / 2.0;
    }

    public Number getY(int s, int i) {
        return getEndY(s, i);
    }

    public double getYValue(int s, int i) {
        return getEndYValue(s, i);
    }

    public int getSeriesCount() {
        return seriesNames.size();
    }

    public Comparable getSeriesKey(int s) {
        return seriesNames.get(s);
    }

    public int indexOf(Comparable s) {
        return seriesNames.indexOf((String)s);
    }

    public void addChangeListener(DatasetChangeListener arg0) {
        logger.log(Level.FINEST, "addChangeListener called");
    }

    public DatasetGroup getGroup() {
        return null;
    }

    public void removeChangeListener(DatasetChangeListener arg0) {
        logger.log(Level.FINEST, "removeChangeListener called");
    }

    public void setGroup(DatasetGroup arg0) {
    } 
}