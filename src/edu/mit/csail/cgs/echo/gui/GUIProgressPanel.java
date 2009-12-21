package edu.mit.csail.cgs.echo.gui;

import java.util.*;

import java.awt.*;

import javax.swing.*;

import edu.mit.csail.cgs.utils.Listener;

public class GUIProgressPanel extends JPanel implements Listener<ChangedEvent> {

	private ProgressPanelInterface.ProcessorPanelImpl progress;
	
	public GUIProgressPanel() { 
		super();
		progress = new ProgressPanelInterface.ProcessorPanelImpl();
		progress.addEventListener(this);
	}
	
	public ProgressPanelInterface getInterface() { return progress; }

	public void eventRegistered(ChangedEvent e) {
		repaint();
	}

	protected void paintComponent(Graphics g) { 
		super.paintComponent(g);
		int w = getWidth(), h = getHeight();
		g.setColor(Color.white);
		g.fillRect(0, 0, w, h);
		
		g.setColor(Color.black);
		HashSet<String> ids = new HashSet<String>(progress.getCurrentIDs());
		TreeSet<String> total = new TreeSet<String>(progress.getTotalIDs());
		
		Graphics2D g2 = (Graphics2D)g;
		FontMetrics fm = g2.getFontMetrics();
		
		int y = 0;
		int x = 1;
		
		int height = 1;
		for(String id : total) { 
			Rectangle r = fm.getStringBounds(id, g).getBounds();
			height = Math.max(height, r.height);
		}

		for(String id : total) { 
			y += (height+3);
			
			g.setColor(ids.contains(id) ? Color.green : Color.red);
			g.fillOval(x, y-height, height, height);
			
			g.setColor(Color.black);
			int c = progress.getCount(id);
			g.drawString("(" + c + ") " + id + " \"" + 
					progress.getDescription(id) + "\"", 
					x+height+3, y);
		}
	}
	
	public Dimension getPreferredSize() { 
		return new Dimension(150, 300);
	}
	
	public static class Frame extends JFrame {
		
		private GUIProgressPanel panel;
		
		public Frame() { 
			super("Echo Calculations");
			panel = new GUIProgressPanel();
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(panel, BorderLayout.CENTER);
			
			setVisible(true);
			pack();
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		}
		
		public GUIProgressPanel getProgressPanel() { return panel; }
	}
}
