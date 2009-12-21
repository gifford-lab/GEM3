/*
 * Created on Aug 16, 2007
 */
package edu.mit.csail.cgs.viz.charting;

import java.util.Collection;
import java.util.Vector;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.event.ChartChangeEvent;
import org.jfree.chart.event.ChartChangeListener;

import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;

/**
 * @author tdanford
 *
 * Just wraps a JFreeChart object up into a Paintable implementation, for integration 
 * in our own modular visualization framework.
 */
public class JFreeChartPaintableAdapter extends AbstractPaintable implements ChartChangeListener {
    
    private Vector<JFreeChart> charts;
    
    public JFreeChartPaintableAdapter() { 
        charts = new Vector<JFreeChart>();
    }
    
    public JFreeChartPaintableAdapter(JFreeChart c) { 
        charts = new Vector<JFreeChart>();
        charts.add(c);
        c.addChangeListener(this);
    }
    
    public JFreeChartPaintableAdapter(Collection<JFreeChart> cs) { 
        charts = new Vector<JFreeChart>();
        for(JFreeChart c : cs) { 
            charts.add(c);
            c.addChangeListener(this);            
        }
    }

    public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
        int w = x2 - x1, h = y2 - y1;
        int ch = (int)Math.floor((double)h / (double)Math.max(1, charts.size()));
        Graphics2D g2 = (Graphics2D)g;
        int cy1 = y1;

        for(JFreeChart chart : charts) { 
            Rectangle2D rect = new Rectangle2D.Double((double)x1, (double)cy1, (double)w, (double)ch);
            chart.draw(g2, rect);
            cy1 += ch;
        }
    }
    
    public void addChart(JFreeChart chart) { 
        charts.add(chart);
        chart.addChangeListener(this);
        dispatchChangedEvent();        
    }

    public void chartChanged(ChartChangeEvent cce) {
        dispatchChangedEvent();
    }
}
