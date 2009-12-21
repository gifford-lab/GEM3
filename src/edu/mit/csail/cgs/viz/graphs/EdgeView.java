/**
 * 
 */
package edu.mit.csail.cgs.viz.graphs;

import java.util.*;
import java.awt.*;
import java.awt.geom.*;

/**
 * @author Timothy Danford
 */
public class EdgeView extends ObjectView {
	
	private GraphView graph;
	private NodeView start, finish;
	
	public EdgeView(GraphView g, NodeView s, NodeView f) {
		super();
		graph = g;
		start = s; finish = f;
		if(start.equals(finish)){setOption("self", new Boolean(true));}
		else{setOption("self", new Boolean(false)); setOption("selfAngle", new Double(Math.PI));}		
        setOption("directed", "false");
        setOption("arrowSize", 10);
	}

	public EdgeView(ObjectView defs, GraphView g, NodeView s, NodeView f) {
		super(defs);
		graph = g;
		start = s; finish = f;
		if(start.equals(finish)){setOption("self", new Boolean(true));}
		else{setOption("self", new Boolean(false));setOption("selfAngle", new Double(Math.PI));}
        setOption("directed", "false");
        setOption("arrowSize", 10);
	}

	public GraphView getGraph() { return graph; }
	public NodeView getStart() { return start; }
	public NodeView getFinish() { return finish; }
    
    public boolean hasDynamicAttributes() { 
        return containsOption("dynamicAttrs");
    }
    
    public void setDynamicAttributes() { 
        setOption("dynamicAttrs", "true");
        
        if(isDirected()) { 
        	int sx = start.getX(), ex = finish.getX();
            int sy = start.getY(), ey = finish.getY();
            
            double dx = (double)(sx - ex), dy = (double)(sy - ey);
            double length = Math.sqrt((dx * dx) + (dy * dy));
            
            double rot = 0.0;
            if(dy >= 0.0) { 
                rot = Math.acos(dx/length);
            } else { 
                rot = -Math.acos(dx/length);
            }
            
            setOption("arrowRotation", rot);
            
            int ewidth = finish.containsOption("width") ? (Integer)finish.getOption("width") : 10;
            int arrowSize = 10;
            
            double rotdiff = Math.PI / 20.0;
            
            Point p1 = rotate(new Point(ewidth, 0), rot);
            Point p2 = rotate(new Point(ewidth+arrowSize,0), rot+rotdiff);
            Point p3 = rotate(new Point(ewidth+arrowSize,0), rot-rotdiff);
            
            p1.translate(ex, ey);
            p2.translate(ex, ey);
            p3.translate(ex, ey);
            
            setOption("arrowPoint1", p1);
            setOption("arrowPoint2", p2);
            setOption("arrowPoint3", p3);
        }
    }
    
    private Point rotate(Point p, double theta) { 
        double cosTheta = Math.cos(theta), sinTheta = Math.sin(theta);
        double xp = (double)p.x * cosTheta - (double)p.y * sinTheta;
        double yp = (double)p.x * sinTheta + (double)p.y * cosTheta;
        return new Point((int)Math.round(xp), (int)Math.round(yp));
    }
    
    public void clearDynamicAttributes() { 
        clearOption("dynamicAttrs");
        
        if(isDirected()) { 
            clearOption("arrowRotation");
            clearOption("arrowPoint1");
            clearOption("arrowPoint2");
            clearOption("arrowPoint3");
        }
    }
    
    public void setDirected(boolean directed) { 
        setOption("directed", String.valueOf(directed));
    }
    
    public boolean isDirected() { 
        String value = getOption("directed").toString();
        if(!value.equals("true") && !value.equals("false")) { 
            throw new IllegalStateException(String.format("option 'directed' has illegal value '%s'",
                    value));
        }
        return value.equals("true");
    }
    
    public void twistSelfEdge(double angle){
    	double currAngle = (Double)getOption("selfAngle");
    	setOption("selfAngle", currAngle+angle);
    }
    
    public void paintView(Graphics2D g2) { 
		if(containsOption("self") && (Boolean)getOption("self")==true) {
			Double angle = (Double)getOption("selfAngle");
			if(angle == null) { angle = Math.PI / 2.0; }
			paintSelfEdge(g2, angle);
		}else if(containsOption("curved")&& (Boolean)getOption("curved")==true) { 
			paintCurvedEdge(g2);
		} else { 
			paintEdge(g2);	
		}				
    }
    
    public void paintEdge(Graphics2D g2) { 
        NodeView first = getStart(), last = getFinish();
        int sr = first.getWidth()/2;
    	int er = last.getWidth()/2;
        int sx = first.getX(), sy = first.getY();
        int ex = last.getX(), ey = last.getY(); 
        
        double omega = Math.atan2((sy-ey),(sx-ex));
        int textOffset=10;

        Stroke oldStroke = g2.getStroke();
        
        if(containsOption("edgeWidth")) { 
            int ewidth = (Integer)getOption("edgeWidth");
            g2.setStroke(new BasicStroke((float)ewidth));
            sx = first.getX() -(new Double((ewidth/2)*Math.cos(omega)).intValue());
            ex = last.getX() +(new Double((ewidth/2)*Math.cos(omega)).intValue());
            sy = first.getY() -(new Double((ewidth/2)*Math.sin(omega)).intValue());
            ey = last.getY() +(new Double((ewidth/2)*Math.sin(omega)).intValue());
        }
        
        g2.setColor(getColor());
        g2.drawLine(sx, sy, ex, ey);
        
        int[] xs = new int[3], ys = new int[3];
        
        if(isDirected() && (ex != sx || ey != sy)) { 
            if(!hasDynamicAttributes()) { 
                setDynamicAttributes();
            }
            
            //Tim's way of drawing arrows:
            /*Point p1 = (Point)ev.getOption("arrowPoint1");
            Point p2 = (Point)ev.getOption("arrowPoint2");
            Point p3 = (Point)ev.getOption("arrowPoint3");
            xs[0] = p1.x; xs[1] = p2.x; xs[2] = p3.x; 
            ys[0] = p1.y; ys[1] = p2.y; ys[2] = p3.y; 
            g2.fillPolygon(xs, ys, 3);
            if(ev.containsOption("name")){
				String edName = (String)ev.getOption("name");
				g2.setColor(Color.black);
				g2.drawString(edName, xs[2]+textOffset, ys[2]);
				g2.setColor(ev.getColor());
			}*/
            
            //Shaun's way of drawing arrows:
            int c0 = (int)Math.round(((ex+sx)/2));
			int c1 = (int)Math.round(((ey+sy)/2));
			
			GraphPaintable.drawArrow(g2, (Integer)getOption("arrowSize"), 
					new Point(c0, c1), omega+Math.PI);
			
			if(containsOption("name")){
				String edName = (String)getOption("name");
				g2.setColor(Color.black);
				g2.drawString(edName, c0+textOffset, c1);
				g2.setColor(getColor());
			}
        }
        
        g2.setStroke(oldStroke);
    }
    
    public void paintCurvedEdge(Graphics2D g2) { 
        NodeView first = getStart(), last = getFinish();
        int sx = first.getX(), sy = first.getY();
        int ex = last.getX(), ey = last.getY();
        double omega = Math.atan2((sy-ey),(sx-ex));
        int ewidth=10;
        Stroke oldStroke = g2.getStroke();
        
        if(containsOption("edgeWidth")) { 
            ewidth = (Integer)getOption("edgeWidth");
            g2.setStroke(new BasicStroke((float)ewidth));
        }
        
        g2.setColor(getColor());
        
        int curveBend = options.containsKey("curve-bend") ? 
        		(Integer)options.get("curve-bend") : 250;
        
        AffineTransform currTrans = g2.getTransform();
        Point po1 = new Point(sx, sy);
        Point po2 = new Point(ex, ey);
        double vx = 0;
        double vy = curveBend; //Controls how "bent" the curve is
        int textOffset=10; 
        
        g2.translate(sx+(ex-sx)/2, sy+(ey-sy)/2);
        Point nP = rotate(new Point(new Double(vx).intValue(),new Double(vy).intValue()), omega);
       
        QuadCurve2D.Double curve =
        	new QuadCurve2D.Double(-(ex-sx)/2, -(ey-sy)/2, new Double(nP.getX()).intValue(), new Double(nP.getY()).intValue(), (ex-sx)/2, (ey-sy)/2);
        g2.draw(curve);
        if(isDirected()){
			PathIterator itr = curve.getPathIterator(null);
			double[] coords = new double[6];
			double[] prev = new double[6];
			do { 
				int type = itr.currentSegment(coords);
				switch(type) { 
				case PathIterator.SEG_MOVETO:
					break;
					
				case PathIterator.SEG_LINETO:
					break;
					
				case PathIterator.SEG_CUBICTO:
         
				case PathIterator.SEG_QUADTO:
					int c0 = (int)Math.round((coords[0]+coords[4])/2);
					int c1 = (int)Math.round((coords[1]+coords[5])/2);
					
					GraphPaintable.drawArrow(g2, (Integer)getOption("arrowSize"), 
							new Point(c0, c1), omega+Math.PI);
					
					if(containsOption("name")){
						String edName = (String)getOption("name");
						g2.setColor(Color.black);
						g2.drawString(edName, c0+textOffset, c1);
						g2.setColor(getColor());
					}
					break;
					
				case PathIterator.SEG_CLOSE:
					break;
				}
				
				for(int i = 0; i < coords.length; i++) { 
					prev[i] = coords[i];
				}
				itr.next();
			} while(!itr.isDone());
        }
		g2.setTransform(currTrans);
		g2.setStroke(oldStroke);            
    }
    
    public void paintSelfEdge(Graphics2D g2, double angle) { 
    	double offset = 200;
        double offsetAngle = Math.toRadians(50);
        int textOffset=10;
        
        double omega = angle;
        NodeView first = getStart(), last = getFinish();
        int sx = first.getX(), sy = first.getY();
        int ex = last.getX(), ey = last.getY();

        Stroke oldStroke = g2.getStroke();
        
        if(containsOption("edgeWidth")) { 
            int ewidth = (Integer)getOption("edgeWidth");
            g2.setStroke(new BasicStroke((float)ewidth));
        }            
        g2.setColor(getColor());
        
        AffineTransform currTrans = g2.getTransform();
        
        g2.translate(sx, sy);
        Point cp1 = rotate(new Point(0,new Double(offset).intValue()), omega-offsetAngle);
        Point cp2 = rotate(new Point(0,new Double(offset).intValue()), omega+offsetAngle);
        
        CubicCurve2D.Double curve =
        	new CubicCurve2D.Double(0, 0, new Double(cp1.getX()).intValue(), new Double(cp1.getY()).intValue(), new Double(cp2.getX()).intValue(), new Double(cp2.getY()).intValue(), 0, 0);
        
        g2.draw(curve);
        if(isDirected()){
        	PathIterator itr = curve.getPathIterator(null);
        
			double[] coords = new double[6];
			double[] prev = new double[6];
			do { 
				int type = itr.currentSegment(coords);
				switch(type) { 
				case PathIterator.SEG_MOVETO:
					break;
					
				case PathIterator.SEG_LINETO:
					break;
					
				case PathIterator.SEG_QUADTO:         
				case PathIterator.SEG_CUBICTO:
					int c0 = (int)Math.round((coords[0]+coords[2])/2);
					int c1 = (int)Math.round((coords[1]+coords[3])/2);
					int c2 = (int)Math.round(3*c0/4);
					int c3 = (int)Math.round(3*c1/4);
					
					GraphPaintable.drawArrow(g2, (Integer)getOption("arrowSize"), 
							new Point(c2, c3), omega+Math.PI);
					
					if(containsOption("name")){
						String edName = (String)getOption("name");
						g2.setColor(Color.black);
						g2.drawString(edName, c2+textOffset, c3);
						g2.setColor(getColor());
					}
					break;
					
				case PathIterator.SEG_CLOSE:
					break;
				}
				
				for(int i = 0; i < coords.length; i++) { 
					prev[i] = coords[i];
				}
				itr.next();
			} while(!itr.isDone());
        }    
		g2.setTransform(currTrans);g2.setColor(Color.red); g2.fillOval(sx, sy, 4, 4);
		g2.setStroke(oldStroke);            
    }
}
