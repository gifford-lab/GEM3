/**
 * 
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.*;
import javax.swing.*;
import java.awt.*;


/**
 * @author Timothy Danford
 */
public class MessageFrame extends JFrame {

	public MessageFrame(String m) {
		super("Message");
		JLabel lbl = new JLabel(m);
		Container c= (Container)getContentPane();
		c.setLayout(new BorderLayout());
		c.add(lbl, BorderLayout.NORTH);
		
		setLocation(50, 50);
		setVisible(true);
		pack();
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}

}
