/**
 * 
 */
package edu.mit.csail.cgs.viz.graphs;

import java.awt.*;
import java.util.*;

/**
 * @author Timothy Danford
 */
public class NodeView extends ObjectView {
	
	public static enum Shape { CIRCLE, SQUARE, HEXAGON };
	
	private GraphView graph;
    private Vector<EdgeView> outgoingEdges;
    private Vector<EdgeView> incomingEdges;

	public NodeView(GraphView g) {
		super();
		graph = g;
        outgoingEdges = new Vector<EdgeView>();
        incomingEdges = new Vector<EdgeView>();
	}
	
	public NodeView(ObjectView defs, GraphView g) {
		super(defs);
		graph = g;
		outgoingEdges = new Vector<EdgeView>();
        incomingEdges = new Vector<EdgeView>();
        options.put("width", 10);
	}
	
	public boolean containsPoint(Point p) { 
		int hh = getHeight() / 2, hw = getWidth()/2;
		return Math.abs(p.x - getX()) <= hw && Math.abs(p.y - getY()) <= hh;
	}
	
	public void paintView(Graphics2D g2) { 
		
		int x = getX(), y = getY();
		int width = getWidth();
		int height = getHeight();
		
		Shape shape = getShape();
		
		int h = width, w = height;
		
		if(containsName()) { 
			Font font = getFont();
			Font oldFont = g2.getFont();
			g2.setFont(font);
			
			String n = getName();
			FontMetrics fm = g2.getFontMetrics();
			
			int padding = 4;
			h = fm.getAscent() + fm.getDescent() + padding;
			w = fm.charsWidth(n.toCharArray(), 0, n.length()) + padding;
			
			w *= 1.1;
			h *= 1.3;
			
			setWidth(w);
			setHeight(h);
			
			g2.setFont(font);
		}
		
		if(shape.equals(Shape.CIRCLE)) { 
			g2.setColor(Color.white);
			//g2.fillOval(x-w/2, y-h/2, w, h);
			g2.setColor(getColor());
			g2.drawOval(x-w/2, y-h/2, w, h);
		} else if (shape.equals(Shape.SQUARE)) { 
			g2.setColor(Color.white);
			//g2.fillRect(x-w/2, y-h/2, w, h);
			g2.setColor(getColor());
			g2.drawRect(x-w/2, y-h/2, w, h);
		} else if (shape.equals(Shape.HEXAGON)) { 
			g2.setColor(Color.white);
			g2.fillRect(x-w/2, y-h/2, w, h);
			g2.setColor(getColor());
			g2.drawRect(x-w/2, y-h/2, w, h);
			g2.drawOval(x-w/2, y-h/2, w, h);
		}
	}
	
	public void paintName(Graphics2D g2) {
		if(containsOption("name")) { 
			int x = getX(), y = getY();
			int w = getWidth(), h = getHeight();

			g2.setColor(Color.blue);
			Font font = g2.getFont();
			Font newFont = getFont();
			g2.setFont(newFont);

			FontMetrics fm = g2.getFontMetrics();
			String name = getName();
			int nw = fm.charsWidth(name.toCharArray(), 0, name.length());
			int nh = fm.getAscent() + fm.getDescent();
			
			int xOffset = containsName() ? -nw/2 : getTextOffset();
			//int yOffset = containsName() ? h/2-fm.getDescent()-1 : xOffset;
			
			int textHeight = fm.getAscent();
			int yOffset = containsName() ? h/2 - textHeight/2 : xOffset;
			
			g2.drawString(name, x + xOffset, y + yOffset);
			
			g2.setFont(font);
		}
	}
	
	public Font getFont() {
        if(containsOption("font")) { 
            return (Font)options.get("font");
        } else { 
            return new Font("Arial", Font.BOLD, 16);
        }
	}
    
    public void setFont(Font f) { 
        options.put("font", f);
    }

	public void setShape(Shape s) { 
		options.put("shape", s);
	}
	
	public Shape getShape() { 
		if(containsOption("shape")) { 
			return (Shape)options.get("shape"); 
		} else { 
			return Shape.CIRCLE;
		}
	}
	
	public boolean containsName() { 
		return containsOption("contains-name") && 
			(Boolean)options.get("contains-name") == true;
	}
	
	public void setContainsName(boolean v) { 
		options.put("contains-name", v);
	}
	
	public int getTextOffset() { 
		return containsOption("text-offset") ? (Integer)options.get("text-offset") : 20;
	}
	
	public void setTextOffset(int to) { 
		options.put("text-offset", to);
	}
	
	public int getNumEdges() { return outgoingEdges.size(); }
	public EdgeView getEdge(int i) { return outgoingEdges.get(i); }
	public GraphView getGraph() { return graph; }
	
	public int getWidth() { 
		if(containsOption("width")) { 
			return (Integer)options.get("width"); 
		} else { 
			return 10;
		}
	}
	
	public int getHeight() { 
		if(containsOption("height")) { 
			return (Integer)options.get("height");
		} else { 
			return getWidth();
		}
	}
	
	public void setWidth(int w) { 
		options.put("width", w);
	}
	
	public void setHeight(int h) { 
		options.put("height", h);
	}
	
	public String getTimeMark() { return (String)options.get("timemark"); }
	public void setTimeMark(String n) { options.put("timemark", n); }
	public void addEdge(EdgeView ev) { 
        outgoingEdges.add(ev);
        ev.getFinish().incomingEdges.add(ev);
	}
    
    public int getNumIncomingEdges() { return incomingEdges.size(); }
    public EdgeView getIncomingEdge(int i) { return incomingEdges.get(i); }

}
