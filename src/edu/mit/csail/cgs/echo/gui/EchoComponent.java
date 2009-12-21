/*
 * Created on Feb 19, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.font.*;
import java.awt.geom.*;
import java.util.*;

import edu.mit.csail.cgs.echo.*;

public abstract class EchoComponent implements EchoEdgeConnector {
    
    private static final double cos_pi4 = Math.cos(Math.PI / 4.0);

    private Rectangle rect;
    private String name;
    private boolean selected;
    private EchoGUI guiComponent;
    
    protected Map<String,Object> params;
    protected JPopupMenu popupMenu;
    
    public EchoComponent(String n, Rectangle r, EchoGUI g) {
        name = n;
        rect = r;
        params = new HashMap<String,Object>();
        selected = false;
        guiComponent = g;
    
        popupMenu = new JPopupMenu();

        JMenuItem item;
        popupMenu.add(item = new JMenuItem("Erase Edges"));
        item.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                clearEdges();
                repaintGUI();
            }
        });
        
        popupMenu.add(item = new JMenuItem("Remove"));
        item.addActionListener(new ActionListener() {
        	public void actionPerformed(ActionEvent e) {
        		guiComponent.remove(EchoComponent.this);
        	}
        });
    }
    
    protected void repaintGUI() { guiComponent.repaint(); }
    
    public EchoGUI getGUI() { return guiComponent; }
    
    public void move(int dx, int dy) { 
        rect.setLocation(getX() + dx, getY() + dy);
    }
    
    public String getName() { return name; }
    public void setName(String n) { name = n; }
    
    public boolean isSelected() { return selected; }
    public void setSelected(boolean v) { selected = v; }
    public void toggleSelected() { selected = !selected; }
    
    public void setParam(String k, Object v) { params.put(k, v); }

    public Point getULPoint() { return rect.getLocation(); }
    public Point getLRPoint() { return new Point(getX() + getWidth(), getY() + getHeight()); }
    public boolean containsPoint(Point p) { return rect.contains(p); }
    
    public int getX() { return rect.x; }
    public int getY() { return rect.y; }
    public int getWidth() { return rect.width; }
    public int getHeight() { return rect.height; }
    public int getRadius() { return (int)Math.ceil(cos_pi4 * (double)Math.max(rect.width, rect.height) / 2.0); } 
    
    public Point getCenterPoint() { return new Point(getX() + getWidth()/2, getY() + getHeight()/2); }
    
    public JPopupMenu getPopupMenu() { return popupMenu; }
    
    public abstract void paint(Graphics2D g);
    
    public static Color lighten(Color c) { 
    	int r = c.getRed(), g = c.getGreen(), b = c.getBlue();
    	int sum = Math.max(r + g + b, 1);
    	double rf = (double)r / (double)sum;
    	double gf = (double)g / (double)sum;
    	double bf = (double)b / (double)sum;
    	
    	sum = (int)Math.ceil(0.5 * sum);
    	r = Math.min(255, (int)Math.round(rf * sum));
    	g = Math.min(255, (int)Math.round(gf * sum));
    	b = Math.min(255, (int)Math.round(bf * sum));
    	
    	return new Color(r, g, b, c.getAlpha());
    }
}
