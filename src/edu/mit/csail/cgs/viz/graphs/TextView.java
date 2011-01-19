package edu.mit.csail.cgs.viz.graphs;

import java.awt.*;
import java.util.Vector;

import edu.mit.csail.cgs.viz.utils.TextLayout;

public class TextView extends NodeView {
	
	private GraphView graph;
    
	public TextView(GraphView g) {
		super(g);
		graph = g;
	}
	
	public TextView(ObjectView defs, GraphView g) {
		super(defs, g);
		graph = g;
	}
	
	public TextView(ObjectView defs, GraphView g, String para) {
		super(defs, g);
		graph = g;
        options.put("paragraph", para);
	}
	
	public void paintName(Graphics2D g2) {
		int x = getX(), y = getY(); 
		if(containsOption("name")) { 
			g2.setColor(Color.black);
			Font font = g2.getFont();
			Font newFont = new Font(font.getName(), Font.BOLD, 16);
			g2.setFont(newFont);
			g2.drawString(getName(), x, y);
		}
	}
	
	public void paintView(Graphics2D g2) { 
		int x = getX(), y = getY(); 
		if(containsOption("paragraph")) { 
			Font oldFont = g2.getFont();
			Font f = containsOption("font") ? 
					(Font)options.get("font") : 
					new Font("Arial", Font.PLAIN, 12);
			g2.setFont(f);
					
			TextLayout layout = new TextLayout();
			int lineLength = containsOption("line-length") ? 
					(Integer)options.get("line-length") : 40;
					
			String paragraph = (String)options.get("paragraph");
			Vector<String> lines = layout.paragraphLayout(paragraph, lineLength);
			
			FontMetrics fm = g2.getFontMetrics();
			int lineHeight = fm.getAscent() + fm.getDescent();

			int lineWidth = 1;
			for(String line : lines) { 
				int lw = fm.charsWidth(line.toCharArray(), 0, line.length());
				lineWidth = Math.max(lw, lineWidth);
			}

			int totalWidth = lineWidth;
			int totalHeight = (lineHeight+1)*lines.size();
			setWidth(totalWidth);
			setHeight(totalHeight);
			
			int sx = x - totalWidth/2;
			int sy = y - totalHeight/2;
			
			g2.setColor(Color.white);
			g2.fillRect(sx, sy, totalWidth, totalHeight);
			g2.setColor(Color.black);
			//g2.drawRect(sx, sy, totalWidth, totalHeight);
			
			int lx = sx, ly = sy + lineHeight;
			for(int i = 0; i < lines.size(); i++) { 
				g2.drawString(lines.get(i), lx, ly);
				ly += lineHeight + 1;
			}
			
			g2.setFont(oldFont);
		}
	}
}
