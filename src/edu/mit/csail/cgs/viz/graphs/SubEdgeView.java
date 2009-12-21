/**
 * 
 */
package edu.mit.csail.cgs.viz.graphs;

import java.util.*;
import java.awt.*;
import java.awt.geom.*;

/**
 * @author Shaun Mahony
 */
public class SubEdgeView extends EdgeView {
	
	private int totalWidth;
	private int numSubEdges;
	private int subEdgeWidths[];
	private Color subEdgeColors[];
	
	public SubEdgeView(GraphView g, NodeView s, NodeView f) {
		super(g,s,f);
		numSubEdges=0;
		subEdgeWidths = new int [50];
		subEdgeColors = new Color [50];
	}

	public SubEdgeView(ObjectView defs, GraphView g, NodeView s, NodeView f) {
		super(defs,g,s,f);
		numSubEdges=0;
		subEdgeWidths = new int [50];
		subEdgeColors = new Color [50];
	}

	public int getTotalWidth(){return totalWidth;}
	public int getNumSubEdges(){return numSubEdges;}
	public int getSubWidth(int x){if(x>=0 && x<numSubEdges){return subEdgeWidths[x];}else{return -1;}}
	public Color getSubColor(int x){if(x>=0 && x<numSubEdges){return subEdgeColors[x];}else{return Color.black;}}
	
	public void setTotalWidth(int w){totalWidth=w;}
	public void addSubEdge(int w, Color c){
		subEdgeWidths[numSubEdges]=w;
		subEdgeColors[numSubEdges]=c;
		numSubEdges++;		
	}
	
	public void reset(){
		numSubEdges=0;
	}
	
	public void paintView(Graphics2D g2) { 
		if(containsOption("self") && (Boolean)getOption("self")==true) {
			double angle = (Double)getOption("selfAngle");
			paintSelfSubEdges(g2, angle);
			
		}else if(containsOption("curved")&& (Boolean)getOption("curved")==true) { 
			paintCurvedSubEdges(g2);
			
		} else { 
			paintSubEdges(g2);	
		}		
	}
	
    public void paintSubEdges(Graphics2D g2) { 
        NodeView first = getStart(), last = getFinish();
        int sr = first.getWidth()/2;
    	int er = last.getWidth()/2; 
        int sx = first.getX(), sy = first.getY();
        int ex = last.getX(), ey = last.getY();
        double omega = Math.atan2((sy-ey),(sx-ex));
        Stroke oldStroke = g2.getStroke();
        
        int totalWidth = getTotalWidth();
        int cumul=0;
        for(int i=0; i < getNumSubEdges(); i++){
        	int ewidth = getSubWidth(i);
        	g2.setStroke(new BasicStroke((float)ewidth));
            int offset = (totalWidth/2)-(cumul)-(ewidth/2)+1;
            
            sx = first.getX() -(new Double((ewidth/2)*Math.cos(omega)).intValue());
            ex = last.getX() +(new Double((ewidth/2)*Math.cos(omega)).intValue());
            sy = first.getY() -(new Double((ewidth/2)*Math.sin(omega)).intValue());
            ey = last.getY() +(new Double((ewidth/2)*Math.sin(omega)).intValue());
 
        	sx = sx + (new Double(offset*Math.sin(omega)).intValue());
            ex = ex + (new Double(offset*Math.sin(omega)).intValue());
            sy = sy - (new Double(offset*Math.cos(omega)).intValue());
            ey = ey - (new Double(offset*Math.cos(omega)).intValue());
            g2.setColor(getSubColor(i));
            g2.drawLine(sx, sy, ex, ey);
            cumul+=ewidth;
        }
        
        //Sub edges cannot be directed
        
        g2.setStroke(oldStroke);
    }
    
    public void paintCurvedSubEdges(Graphics2D g2) { 
        NodeView first = getStart(), last = getFinish();
        int sx = first.getX(), sy = first.getY();
        int ex = last.getX(), ey = last.getY();
        double omega = Math.atan2((sy-ey),(sx-ex));
        
        Stroke oldStroke = g2.getStroke();
        AffineTransform currTrans = g2.getTransform();
        
        int totalWidth = getTotalWidth();
        int cumul=0;
        for(int i=0; i < getNumSubEdges(); i++){
        	sx = first.getX(); sy = first.getY();
            ex = last.getX(); ey = last.getY();
            
        	int ewidth = getSubWidth(i);
        	g2.setStroke(new BasicStroke((float)ewidth));
        	int offset = (totalWidth/2)-(cumul)-(ewidth/2);
        	
        	double vx = 0;
            double vy = 250; //Controls how "bent" the curve is
            Point oldCent = new Point(sx+(ex-sx)/2, sy+(ey-sy)/2);
            g2.translate(sx+(ex-sx)/2, sy+(ey-sy)/2);
            
            Point nP = GraphPaintable.rotate(
            		new Point(new Double(vx).intValue(),new Double(vy).intValue()), omega);
            
            double sPhi = Math.atan2(new Double(nP.getY()).intValue()+(ey-sy)/2, new Double(nP.getX()).intValue()+(ex-sx)/2)-Math.toRadians(90);
            double ePhi = Math.atan2(new Double(nP.getY()).intValue()-(ey-sy)/2, new Double(nP.getX()).intValue()-(ex-sx)/2)+Math.toRadians(90);
            g2.setTransform(currTrans);
            sx = sx +(new Double(offset*Math.cos(sPhi)).intValue());
            ex = ex +(new Double(offset*Math.cos(ePhi)).intValue());
            sy = sy +(new Double(offset*Math.sin(sPhi)).intValue());
            ey = ey +(new Double(offset*Math.sin(ePhi)).intValue());
            Point newCent = new Point(sx+(ex-sx)/2, sy+(ey-sy)/2);
            double old2newDist = oldCent.distance(newCent);
            g2.translate(sx+(ex-sx)/2, sy+(ey-sy)/2);
            double vyb = ((125-old2newDist+totalWidth/2-(ewidth/2)-(cumul))*2)+1;
            
            Point nPb = GraphPaintable.rotate(
            		new Point(new Double(vx).intValue(),new Double(vyb).intValue()), omega);
            
            QuadCurve2D.Double curve =
            	new QuadCurve2D.Double(-(ex-sx)/2, -(ey-sy)/2, new Double(nPb.getX()).intValue(), new Double(nPb.getY()).intValue(), (ex-sx)/2, (ey-sy)/2);
            
            g2.setColor(getSubColor(i));
            cumul+=ewidth;
            g2.draw(curve);
            g2.setTransform(currTrans);	            
        }	        
                    
        g2.setTransform(currTrans);
		g2.setStroke(oldStroke);            
    }
    
    public void paintSelfSubEdges(Graphics2D g2, double angle) { 
    	double offset = 200;
        double offsetAngle = Math.toRadians(50);
        double omega = angle;
        NodeView first = getStart(), last = getFinish();
        int sx = first.getX(), sy = first.getY();
        int ex = first.getX(), ey = first.getY();
        
        Stroke oldStroke = g2.getStroke();
        
        if(containsOption("edgeWidth")) { 
            int ewidth = (Integer)getOption("edgeWidth");
            g2.setStroke(new BasicStroke((float)ewidth));
        }
        int totalWidth = getTotalWidth();
        int cumul=0;
        for(int i=0; i < getNumSubEdges(); i++){
        	sx = first.getX(); sy = first.getY();
            ex = last.getX(); ey = last.getY();
                           
        	int ewidth = getSubWidth(i);
        	g2.setStroke(new BasicStroke((float)ewidth));
        	int startOffset = (totalWidth/2)-(cumul)-(ewidth/2);
        	
        	
        	AffineTransform currTrans = g2.getTransform();
        	g2.translate(sx, sy);

        	Point basePoint = GraphPaintable.rotate(new Point(-startOffset, -startOffset), omega+Math.toRadians(45));
            Point cpmain = GraphPaintable.rotate(new Point(0,new Double(offset).intValue()), omega);
            Point cpm1 = GraphPaintable.rotate(new Point(0,new Double(offset).intValue()), omega-offsetAngle);
            Point cpm2 = GraphPaintable.rotate(new Point(0,new Double(offset).intValue()), omega+offsetAngle);
            
            int cpmidx = (int)Math.round((cpm1.getX()+cpm2.getX())/2);
            int cpmidy = (int)Math.round((cpm1.getY()+cpm2.getY())/2);
            int newapexx = (cpmidx*3/4) - new Double((((totalWidth-ewidth)/2)-cumul)*Math.sin(omega)).intValue();
            int newapexy = (cpmidy*3/4) + new Double((((totalWidth-ewidth)/2)-cumul)*Math.cos(omega)).intValue();
            double apexRatio = cpmain.distance(0,0)/(new Point(cpmidx,cpmidy).distance(0,0));
            int newcpmainX =new Double((newapexx*4/3)*apexRatio).intValue();
            int newcpmainY =new Double((newapexy*4/3)*apexRatio).intValue();
            int newDist = new Double(new Point(newcpmainX, newcpmainY).distance(basePoint.getX(), basePoint.getY())).intValue();
          
            Point cp1 = GraphPaintable.rotate(new Point(0, newDist), omega-offsetAngle);
            Point cp2 = GraphPaintable.rotate(new Point(0, newDist), omega+offsetAngle);
            
            CubicCurve2D.Double curve =
            	new CubicCurve2D.Double(new Double(basePoint.getX()).intValue(), new Double(basePoint.getY()).intValue(), new Double(cp1.getX()).intValue(), new Double(cp1.getY()).intValue(), new Double(cp2.getX()).intValue(), new Double(cp2.getY()).intValue(), new Double(basePoint.getX()).intValue(), new Double(basePoint.getY()).intValue());

        	g2.setColor(getSubColor(i));
	            cumul+=ewidth;
	            g2.draw(curve);
	            g2.setTransform(currTrans);	            
	        }	                
		g2.setStroke(oldStroke);            
    }
}
