/*
 * Author: tdanford
 * Date: Aug 19, 2008
 */
package edu.mit.csail.cgs.metagenes.swing;

import java.util.*;
import java.util.regex.*;
import java.io.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import edu.mit.csail.cgs.metagenes.*;
import edu.mit.csail.cgs.viz.paintable.PaintableChangedEvent;
import edu.mit.csail.cgs.viz.paintable.PaintableChangedListener;
import edu.mit.csail.cgs.viz.paintable.PaintableScale;
import edu.mit.csail.cgs.viz.utils.FileChooser;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.ewok.verbs.GenomeExpander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;

/**
 * MetaFrame is a Swing component that provides a default interface to the metagene system.
 * 
 * Once it's created, it provides UI elements that allow the user to load new points into 
 * the given profile, and displays the results.  
 * 
 * It also (through calls to the getHandler() method) provides a programmatic hook for adding
 * new Point objects automatically.  
 *  
 * @author tdanford
 */
public class MetaFrame extends JFrame {
	
	private Genome genome;
	private BinningParameters params;
	private MetaProfile profile=null;
	private MetaProfileHandler handler;
	private MetaUtils utils;
	private PaintableScale peakScale, lineScale;
	private ProfileLinePanel linePanel;
	private JButton stopAdding, clusterButton, optionsButton;
	private ProfilePanel panel;
	private OptionsFrame of;
	
	public MetaFrame(Genome g, BinningParameters bps, PointProfiler pp, boolean normalizedMeta) {
		peakScale = new PaintableScale(0.0, 0.0);
		lineScale = new PaintableScale(0.0, 0.0);
		
		genome = g;
		params = bps;
		handler = new MetaProfileHandler("MetaProfile", params, pp, normalizedMeta);
		profile = handler.getProfile();
		linePanel = new ProfileLinePanel(params, lineScale);
		utils = new MetaUtils(genome);
		
		profile.addProfileListener(linePanel);
		
		panel = new ProfilePanel(profile, peakScale);
		
		stopAdding = new JButton(createStopAddingAction());
		clusterButton = new JButton(createClusterAction());
		optionsButton = new JButton(createOptionsAction());
		
		JPanel buttons = new JPanel();
		buttons.setLayout(new FlowLayout());
		buttons.add(stopAdding);
		buttons.add(clusterButton);
		buttons.add(optionsButton);
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		c.add(panel, BorderLayout.CENTER);
		c.add(new JScrollPane(linePanel, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED), BorderLayout.EAST);
		c.add(buttons, BorderLayout.SOUTH);
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setJMenuBar(createJMenuBar());
		
		of = new OptionsFrame(panel, linePanel);
		of.startup();
	}
	
	public MetaProfileHandler getHandler() { return handler; }
	
	private JMenuBar createJMenuBar() { 
		JMenuBar bar = new JMenuBar();
		JMenu menu; JMenuItem item;
		
		bar.add(menu = new JMenu("File"));
		menu.add(item = new JMenuItem(createLoadPointsFileAction()));
		menu.add(item = new JMenuItem(createLoadTSSsAction()));
		menu.add(new JSeparator());
		menu.add(item = new JMenuItem(createExitAction()));
		
		bar.add(menu = new JMenu("Image"));
		menu.add(item = new JMenuItem(panel.createSaveImageAction()));
		menu.add(item = new JMenuItem(linePanel.createSaveImageAction()));
		return bar;
	}
	
	public Action createClusterAction() { 
		return new AbstractAction("Cluster Profiles") { 
			public void actionPerformed(ActionEvent e) { 
				linePanel.cluster();
			}
		};
	}
	
	public Action createStopAddingAction() { 
		return new AbstractAction("Stop Adding...") { 
			public void actionPerformed(ActionEvent e) { 
				handler.stopAllAddingThreads();
			}
		};
	}

	public Action createOptionsAction() { 
		return new AbstractAction("Options...") { 
			public void actionPerformed(ActionEvent e) { 
				of.setVisible(true);
			}
		};
	}
	
	private Action createLoadPointsFileAction() { 
		return new AbstractAction("Load Points...") { 
			public void actionPerformed(ActionEvent e) { 
				FileChooser chooser = new FileChooser(MetaFrame.this);
				File f = chooser.choose();
				if(f != null) { 
					try {
						Vector<Point> points = utils.loadPoints(f);
						handler.addPoints(points);
					} catch (IOException e1) {
						e1.printStackTrace();
					}
				}
			}
		};
	}

	private Action createLoadTSSsAction() { 
		return new AbstractAction("Load TSSs...") { 
			public void actionPerformed(ActionEvent e) { 
				Iterator<Point> points = utils.loadTSSs();
				handler.addPoints(points);
			}
		};
	}
	
	private Action createExitAction() { 
		return new AbstractAction("Exit") { 
			public void actionPerformed(ActionEvent e) {
				MetaFrame.this.dispose();
			}
		};
	}
	public MetaUtils getUtils(){return utils;}
	public void setColor(Color c){
		panel.updateColor(c);
		linePanel.updateColor(c);
	}
	public void setLineMin(double m){linePanel.setMinColorVal(m);}
	public void setLineMax(double m){linePanel.setMaxColorVal(m);}
	public void setLineThick(int t){linePanel.updateLineWeight(t);}
	public void setLinePanelColorQuanta(double [] q){
		linePanel.setLineColorQuanta(q);
	}
	public void startup() {
		// We want to call setVisible() and *then* pack() -- but once setVisible() has 
		// been called on a Swing component, you shouldn't call any other methods of that 
		// component except *from the Swing thread*.  Therefore, this hack.  
		Runnable r = new Runnable() {
			public void run() { 
				MetaFrame.this.setVisible(true);
				MetaFrame.this.pack();
			}
		};
		SwingUtilities.invokeLater(r);
	}
	
	
}

