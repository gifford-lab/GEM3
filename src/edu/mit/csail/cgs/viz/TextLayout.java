/*
 * Author: tdanford
 * Date: Sep 28, 2008
 */
package edu.mit.csail.cgs.viz;

import java.awt.*;
import java.util.*;
import java.util.regex.*;

/**
 * Takes a string, and lays it out in a graphical format, in lines.  We try to pay as much attention 
 * to linebreaks and whitespace as we can, and this class currently performs line-wrapping, although
 * maybe that should be an option to turn on and off.  
 * 
 * @author tdanford
 */
public class TextLayout {
	
	private Font font;
	
	public TextLayout(int size) { 
		font = new Font("Courier", Font.BOLD, size);
	}

	public void paintText(String txt, Graphics2D g2, int x1, int y1, int x2, int y2) { 
		Font oldFont = g2.getFont();
		g2.setFont(font);
		
		FontMetrics fm = g2.getFontMetrics();
		
		String[] lines = txt.split("\n");
		
		int lineHeight = fm.getAscent() + fm.getDescent();
		int y = y1+lineHeight;
		for(int i = 0; i < lines.length; i++, y += lineHeight) {
			y = paintLine(lines[i], g2, x1, x2, y, lineHeight);
		}
		
		g2.setFont(oldFont);
	}
	
	private static Pattern word = Pattern.compile("^([^\\s]*)(\\s*)(.*)$");
	
	private int paintLine(String line, Graphics2D g2, int x1, int x2, int y, int lineHeight) {
		Point p = new Point(x1, y);
		FontMetrics fm = g2.getFontMetrics();
		
		while(line.length() > 0) { 
			Matcher m = word.matcher(line);
			if(m.matches()) { 
				String chars = m.group(1);
				String sp = m.group(2);
				line = m.group(3);
				
				int wordWidth = fm.charsWidth(chars.toCharArray(), 0, chars.length());
				int spaceWidth = fm.charsWidth(sp.toCharArray(), 0, sp.length());
				
				if(p.x+wordWidth > x2) {
					y += lineHeight;
					p = new Point(x1, y);
				}
				
				g2.drawString(chars, p.x, p.y);
				p = new Point(p.x+wordWidth, p.y);
				
				if(p.x+spaceWidth > x2) { 
					y += lineHeight;
					p = new Point(x1, y);
				} else { 
					g2.drawString(sp, p.x, p.y);
					p = new Point(p.x + spaceWidth, p.y);
				}
			} else { 
				return y;
			}
		}
		
		return y;
	}
}
