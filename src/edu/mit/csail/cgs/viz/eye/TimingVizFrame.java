/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package edu.mit.csail.cgs.viz.eye;

import java.io.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;

import edu.mit.csail.cgs.utils.models.*;
import edu.mit.csail.cgs.utils.models.Timer;
import edu.mit.csail.cgs.utils.models.ModelInput.LineReader;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.viz.eye.*;

public class TimingVizFrame extends JFrame implements Timer {
	
	public static void main(String[] args) { 
		TimingVizFrame frame = new TimingVizFrame();
		for(int i = 0; i < args.length; i++) { 
			File f = new File(args[i]);
			try {
				Iterator<Timing> timings = new ModelInputIterator<Timing>(
						new ModelInput.LineReader<Timing>(Timing.class,
								new FileInputStream(f)));
				frame.addTimings(timings);
				System.out.println(String.format("-> %s", args[i]));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	private ModelScatter scatter;
	
	public TimingVizFrame() { 
		super("Timing Data");
		scatter = new ModelScatter("size", "seconds");
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		PaintablePanel pp = new PaintablePanel(scatter);
		pp.setPreferredSize(new Dimension(400, 200));

		c.add(pp, BorderLayout.CENTER);
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		showMe();
	}
	
	private void showMe() { 
		SwingUtilities.invokeLater(new Runnable() { 
			public void run() { 
				setVisible(true);
				pack();
			}
		});
	}
	
	public void addTiming(Timing t) { 
		scatter.addModel(t);
	}
	
	public void addTimings(Iterator<Timing> ts) { 
		while(ts.hasNext()) { 
			addTiming(ts.next());
		}
	}
}
