package edu.mit.csail.cgs.viz.scatter;

import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.chart.*;
import org.jfree.chart.annotations.XYLineAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.category.*; 
import org.jfree.chart.renderer.xy.*; 
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class ScatterPanel extends JPanel {

    private JFreeChart chart;
    private Color color; 
    private DefaultXYDataset dataset;

    public ScatterPanel (String plotTitle,Dataset2D data,double samplingFraction) {
    	this(plotTitle, data, samplingFraction, new Color((float).2, (float)0,(float)1.0,(float).3));
    }

    public ScatterPanel (String plotTitle, Dataset2D data, double samplingFraction, Color c) {
        color = c;
    	dataset = new DefaultXYDataset();
        double[][] values = new double[2][(int)(data.getCount() * samplingFraction)];
        int count = 0;
        for (int i = 0; i < data.getCount(); i++) {
            double one = data.getVal(0,i);
            double two = data.getVal(1,i);
            if (Double.isNaN(one) || Double.isNaN(two) || Double.isInfinite(one)|| Double.isInfinite(two)) {
                continue;
            }else{
	            if (count >= values[0].length) {
	                break;
	            }
	            values[0][count] = one;
	            values[1][count] = two;
	            count++;
            }
        }
        if (count < values.length) {
            while (count < values.length && count > 0) {
                //values[0][count] = values[0][count-1];
                //values[1][count] = values[1][count-1];
            	values[0][count] = 0;
                values[1][count] = 0;
                count++;
            }            
        }
        dataset.addSeries(plotTitle,values);   
        
        chart = ChartFactory.createScatterPlot(plotTitle,
                                               data.getLabelOne(),
                                               data.getLabelTwo(),
                                               dataset,PlotOrientation.VERTICAL,
                                               false,false,false);
        chart.setAntiAlias(true);
        ChartPanel panel = new ChartPanel(chart);
        
        XYPlot plot= chart.getXYPlot();
    
        ValueAxis dAxis = plot.getDomainAxis();
		ValueAxis rAxis = plot.getRangeAxis();
		double maxbound = Math.max(dAxis.getUpperBound(), rAxis.getUpperBound());
		double minbound = Math.min(dAxis.getLowerBound(), rAxis.getLowerBound());
		dAxis.setUpperBound(maxbound);
		rAxis.setUpperBound(maxbound);
		dAxis.setLowerBound(minbound);
		rAxis.setLowerBound(minbound);
		plot.addAnnotation(new XYLineAnnotation(minbound, minbound, maxbound, maxbound));
		plot.addAnnotation(new XYLineAnnotation(0, minbound, 0, maxbound));
		plot.addAnnotation(new XYLineAnnotation(minbound, 0, maxbound, 0));
		
		XYItemRenderer renderer = ((XYPlot)plot).getRenderer();
        renderer.setSeriesPaint(0, color);
        renderer.setBaseStroke(new BasicStroke((float).5));
        renderer.setSeriesStroke(0,new BasicStroke((float).5));
        renderer.setSeriesShape(0,new Rectangle(1,1));
        
        setPreferredSize(new Dimension(800,800));
        setLayout(new BorderLayout());
        add(panel,BorderLayout.CENTER);
    }
    public JFreeChart getChart(){return(chart);}
}