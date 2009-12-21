/*
 * Created on Aug 16, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.viz.charting;

import java.sql.SQLException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.mit.csail.cgs.datasets.expression.*;
import org.jfree.chart.*;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.*;

public class ExpressionCharter {

    private Logger logger;
    private ExpressionLoader loader;
    
    public ExpressionCharter(ExpressionLoader l) { 
        loader = l;
        logger = Logger.getLogger("edu.mit.csail.cgs.viz.charting.ExpressionCharter");
    }
    
    public JFreeChart createScatterPlot(XYDataset ds, String x, String y) { 
        return ChartFactory.createScatterPlot("Plot", x, y, ds, PlotOrientation.VERTICAL, false, false, false);
    }
    
    public JFreeChart createHistogram(Experiment e, Map<String,Set<String>> subsets, int bins) throws SQLException { 
        IntervalXYDataset ds = createHistogramDataset(e, subsets, bins);
        return ChartFactory.createHistogram(null, "Expression Value", "Counts", 
                ds, PlotOrientation.VERTICAL, false, false, false);
    }
    
    public IntervalXYDataset createHistogramDataset(Experiment e, Map<String,Set<String>> subsets, int bins) throws SQLException {
        Map<String,Object> stats = loader.loadExperimentStatistics(e);
        double min = (Double)stats.get("min");
        double max = (Double)stats.get("max");
        
        logger.log(Level.FINEST, "MIN: " + min);
        logger.log(Level.FINEST, "MAX: " + max);
        
        HistogramDataset hd = new HistogramDataset(min, max, bins);
        
        Collection<ExprMeasurement> ms = loader.loadExperimentMeasurements(e);
        Map<String,Double> valuemap = new HashMap<String,Double>();
        for(ExprMeasurement em : ms) { valuemap.put(em.getProbe().getName(), em.getValue()); }

        if(subsets != null) { 
            for(String seriesKey : subsets.keySet()) { 
                Set<String> probeNames = subsets.get(seriesKey);
                int series = hd.addSeries(seriesKey);

                for(String pn : probeNames) { 
                    hd.addValue(series, valuemap.get(pn));
                }
            }
        } else { 
            int series = hd.addSeries("total");
            int size = 0;
            for(String vn : valuemap.keySet()) { 
                hd.addValue(series, valuemap.get(vn));
                size += 1;
            }
            logger.log(Level.FINEST, "Added " + size + " total expression values.");
        }
     
        return hd;
    }
}
