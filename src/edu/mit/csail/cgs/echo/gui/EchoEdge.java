/*
 * Created on Feb 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.font.*;
import java.util.*;

public class EchoEdge {
    
    public static double twopi = 2.0 * Math.PI;

    private EchoComponent comp1, comp2;
    private String name;
    
    public EchoEdge(EchoComponent c1, EchoComponent c2) { 
        comp1 = c1; 
        comp2 = c2;
        name = null;
    }
    
    public EchoEdge(EchoComponent c1, EchoComponent c2, String n) { 
        comp1 = c1; 
        comp2 = c2;
        name = n;
    }
    
    public String getName() { return name; }
    
    public void paint(Graphics2D g) { 
        g.setColor(Color.black);
        Point p1 = comp1.getCenterPoint(), p2 = comp2.getCenterPoint();

        Font oldFont = g.getFont();
        Font font = new Font("Times New Roman", Font.PLAIN, 14);
        g.setFont(font);
        
        int max1 = Math.max(comp1.getWidth(), comp1.getHeight());
        int max2 = Math.max(comp2.getWidth(), comp2.getHeight());
        
        int dx = p2.x - p1.x, dy = p1.y - p2.y;
        double theta = 0.0;

        if((dx != 0 || dy != 0)) { 
            int cx = (p1.x + p2.x)/2, cy = (p1.y + p2.y)/2;

            double r = Math.sqrt(dx*dx + dy * dy);
            int r2 = (int)Math.round(r / 2.0);
            
            if(dx != 0) {
                if(dx > 0) {  
                    theta = Math.atan((double)dy / (double)dx);
                    theta = twopi - theta;
                } else { 
                    theta = - Math.atan((double)dy / (double)dx) - Math.PI;
                    if(name != null) { theta += Math.PI; }
                }
            } else { 
                theta = dy > 0 ? Math.PI / 2.0 : 3.0 * Math.PI / 2.0;
                theta += Math.PI;
            }

            g.translate(cx, cy);
            g.rotate(theta);
            
            if(name != null) {
            	g.setColor(Color.black);
                g.drawLine(-r2, 0, r2, 0);
                FontRenderContext frc = g.getFontRenderContext();
                TextLayout layout = new TextLayout(name, g.getFont(), frc);
                
                int rad1 = comp1.getRadius(), rad2 = comp2.getRadius();
                int textArea = r2 * 2 - rad1 - rad2;
                int layoutWidth = (int)Math.ceil(layout.getBounds().getWidth());
                
                if((int)Math.ceil(layoutWidth) <= textArea) { 
                    g.drawString(name, -layoutWidth/2, -1);
                }
            } else {
            	
                int depth = 3;
                int sep = depth*2;
                int offset = r2/2;
                
            	/*
                g.drawLine(-r2, depth, r2, depth);
                g.drawLine(-r2, -depth, r2, -depth);
                
                g.drawLine(-sep + offset, 0, -sep-depth + offset, -depth);
                g.drawLine(-sep + offset, 0, -sep-depth + offset, depth);
                g.drawLine(offset, 0, -depth + offset, -depth);
                g.drawLine(offset, 0, -depth + offset, depth);
                g.drawLine(sep + offset, 0, sep-depth + offset, -depth);
                g.drawLine(sep + offset, 0, sep-depth + offset, depth);
                */
            	
            	Stroke oldStroke = g.getStroke();
            	g.setStroke(new BasicStroke((float)3.0));
            	int rad = comp2.getRadius();
            	depth = comp2.getRadius();
            	
            	g.setColor(Color.cyan);
                int ptx = r2-(int)(rad*1.5);
                g.drawLine(-r2, 0, ptx, 0);
                g.drawLine(ptx, 0, ptx-depth, -depth);
                g.drawLine(ptx, 0, ptx-depth, depth);
            	
            	g.setStroke(oldStroke);
            }
            
            g.rotate(-theta);
            g.translate(-cx, -cy);
        }
        
        g.setFont(oldFont);
    }
}
