/**
 * @author tdanford
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;
import java.io.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.paintable.*;

/**
 * @author Timothy Danford
 */
public class ProgressPanel extends JPanel implements ProgressListener, PaintableChangedListener {
	
	private ProgressPaintable pp;
	
	public ProgressPanel() {
		super();
		setLayout(new FlowLayout());
		pp = new ProgressPaintable();
		add(new PaintablePanel(pp));
		
		pp.addPaintableChangedListener(this);
	}

	public Dimension getPreferredSize() { 
		return new Dimension(300, 100);
	}
	
	public void progressMade(ProgressEvent pe) { 
		pp.progressMade(pe);
	}
	
	public void registerEventMax(String key, int max) { 
		pp.registerEventMax(key, max);
	}

	public void paintableChanged(PaintableChangedEvent pce) {
		repaint();
	}

	public static class Frame extends JFrame implements ProgressListener {
		
		private ProgressPanel pp;
		
		public Frame() { 
			super("Progress");
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(pp = new ProgressPanel());
			
			setLocation(100, 100);
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setVisible(true);
			pack();
		}
		
		public ProgressPanel getPanel() { return pp; }
		
		public void progressMade(ProgressEvent pe) { 
			pp.progressMade(pe);
		}

		public void registerEventMax(String key, int max) {
			pp.registerEventMax(key, max);
		}
	}

}

class ProgressPaintable extends AbstractPaintable implements ProgressListener {
	
	private Color c;
	private Map<String,Integer> maxima;
	private Map<String,Integer> values;
	
	public ProgressPaintable() { 
		c = Color.blue;
		maxima = new LinkedHashMap<String,Integer>();
		values = new LinkedHashMap<String,Integer>();
	}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int rows = maxima.size();
		int w = x2-x1, h = y2-y1;
		int rowHeight = h/(rows*2+1);
		int barWidth = 4*w/5;
		int barOffset = (w-barWidth)/2;
		
		g.setColor(Color.white);
		g.fillRect(x1, y1, w, h);
		
		Vector<String> keys = new Vector<String>(maxima.keySet());
		for(int i = 0; i < keys.size(); i++) { 
			String key = keys.get(i);
			
			int y = y1 + rowHeight*(i*2+1);
			double f = Math.min(1.0, (double)values.get(key) / (double)maxima.get(key));
			int bar = (int)Math.round(f * (double)barWidth);
		
			g.setColor(c);
			g.fillRect(barOffset, y, bar, rowHeight);
			g.setColor(Color.black);
			g.drawRect(barOffset, y, barWidth, rowHeight);
			
			g.drawString(key, barOffset+1, y-1);
		}
	}

	public void progressMade(ProgressEvent e) {
		String key = e.getKey();
		if(maxima.containsKey(key)) { 
			values.put(key, e.getValue());
		}
		dispatchChangedEvent();
	}

	public void registerEventMax(String key, int max) {
		maxima.put(key, max);
		values.put(key, 0);
		dispatchChangedEvent();
	}
}
                              