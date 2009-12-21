package edu.mit.csail.cgs.echo.gui;

import edu.mit.csail.cgs.echo.ProcessCounter;
import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.text.*;

public class CountPanel extends JPanel implements ProcessCounter {
	
	private static NumberFormat nf;
	
	static { 
		nf = NumberFormat.getInstance();
		nf.setGroupingUsed(true);
	}
	
	private int count;
	private JLabel label;
	
	public CountPanel() { 
		super();
		count = 0;
		label = new JLabel(nf.format(count));
        setLayout(new BorderLayout());
        add(label, BorderLayout.NORTH);
	}

	public void reset() {
		count = 0;
		rebuildLabel();
	}

	public void tick() {
		count += 1;
		rebuildLabel();
	}

	private void rebuildLabel() {
		label.setText(nf.format(count));
		validate();
	}
	
	public static class Frame extends JFrame {
		
		private CountPanel cp;
		public Frame(CountPanel p) { 
			super("Count"); 
			cp = p;
			Container c = (Container)getContentPane();
			c.setLayout(new BorderLayout());
			c.add(cp, BorderLayout.NORTH);
			
			setVisible(true);
			pack();
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		}
		
	}
}
