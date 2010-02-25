/*
 * Created on Sep 7, 2006
 */
package edu.mit.csail.cgs.viz.expression;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.text.*;

import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.utils.io.parsing.affyexpr.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class AffyExpressionPainter extends AbstractPaintable {
    
    private Vector<AffyExperiment> expts;
    private Vector<AffyProbe> probes;

    public AffyExpressionPainter(Collection<AffyExperiment> aes) {
        expts = new Vector<AffyExperiment>(aes);
        probes = new Vector<AffyProbe>();
    }
    
    public void addAffyProbe(AffyProbe ap) { 
        if(!probes.contains(ap)) { probes.add(ap); }
    }

    public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
        int w = x2 - x1, h = y2 - y1;
        Graphics2D g2 = (Graphics2D)g;
        g2.setColor(Color.white);
        g2.fillRect(x1, y1, w, h);
        
        
    }
}

class CachedProbe { 
    
    private AffyProbe probe;
    private Vector<Double> values;
    private Color color;
    private Stroke stroke;
    private double min, max;
    
    public CachedProbe(AffyProbe ap, Vector<AffyExperiment> expts) { 
        probe = ap;
        min = max = 0.0;
        color = Color.lightGray;
        stroke = new BasicStroke((float)2.0);
        values = new Vector<Double>();
        for(AffyExperiment e : expts) { 
            AffyMeasurement am = e.getMeasurement(probe);
            if(am != null) { 
                values.add(am.getValue());
                min = Math.min(am.getValue(), min);
                max = Math.max(am.getValue(), max);
            } else { 
                values.add(null);
            }
        }
    }
    
    public double getMax() { 
        return max;
    }
    
    public double getMin() { 
        return min;
    }
    
    public void paintProbe(Graphics2D g2, int x1, int y1, int x2, int y2, double tmin, double tmax) { 
        g2.setStroke(stroke);
        g2.setColor(color);
        
        int ppx = -1, ppy = -1;
        for(int i = 0; i < values.size(); i++) { 
            Double v = values.get(i);
            if(v != null && tmax - tmin > 0.0) { 
                double frac = (v - tmin) / (tmax - tmin);
                int pixY = y2 - (int)Math.round(frac * (y2 - y1));
                int pixX = x1 + (int)Math.floor((double)(i+1) / (double)(values.size() + 2));
                
                if(ppx != -1) { 
                    g2.drawLine(ppx, ppy, pixX, pixY);
                }
                
                ppx = pixX; ppy = pixY;
                
            } else { 
                ppx = ppy = -1;
            }
        }
    }
}