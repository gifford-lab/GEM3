package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import java.io.IOException;
import java.sql.SQLException;

import javax.swing.*;

import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.warpdrive.WarpOptions;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;;

public class WarpOptionsFrame extends JFrame implements ActionListener {

	static Logger logger = Logger.getLogger(WarpOptionsFrame.class);

	private String genomeString; // hack to run local version
	private ArrayList<PainterContainer> pcs;
	private WarpOptionsPane pane;
	private JButton ok, cancel;

	// variables for menus
	private JMenuBar menuBar;
	private JMenu fileMenu;
	private JMenuItem openSessionItem;
	private JMenuItem saveSessionItem;
	private JMenuItem exitItem;
	private JMenu toolsMenu;
	private JMenuItem optionsItem;

	public WarpOptionsFrame() throws NotFoundException {
		super();
		setTitle("Warp Drive");
		pcs = new ArrayList<PainterContainer>();
		pane = new WarpOptionsPane();
		init();
	}

	public WarpOptionsFrame(String species, String genome) throws NotFoundException {
		super();
		setTitle("Warp Drive for " + species + ", " + genome);
		pcs = new ArrayList<PainterContainer>();
		pane = new WarpOptionsPane(species, genome);
		init();
	}

	public WarpOptionsFrame(WarpOptions opts) throws NotFoundException {
		super();
		setTitle("Warp Drive");
		pcs = new ArrayList<PainterContainer>();
		pane = new WarpOptionsPane(opts);
		init();
		genomeString = opts.genomeString;
	}

	private void init() {
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new GridBagLayout());
		Dimension buttonSize = new Dimension(30, 20);
		ok = new JButton("OK");
		cancel = new JButton("Cancel");
		ok.setMaximumSize(buttonSize);
		cancel.setMaximumSize(buttonSize);
		ok.addActionListener(this);
		cancel.addActionListener(this);
		buttonPanel.add(ok);
		buttonPanel.add(cancel);
		Container content = getContentPane();
		content.setLayout(new BorderLayout());
		content.add(buttonPanel, BorderLayout.SOUTH);
		content.add(pane, BorderLayout.CENTER);

		WarpOptions options = pane.parseOptions();
		this.setSize(options.getPreferredWindowWidth(), options.getPreferredWindowHeight());
		if (options.isWindowCentered()) {
			this.setLocationRelativeTo(null);
		} else {
			this.setLocation(options.getPreferredTopLeftX(), options.getPreferredTopLeftY());
		}

		this.createMenu();

		setVisible(true);
	}

	/**
	 * Create a JMenuBar for this GUI
	 */
	private void createMenu() {
		menuBar = new JMenuBar();

		// build the file menu
		fileMenu = new JMenu("File");
		fileMenu.setMnemonic(KeyEvent.VK_F);
		menuBar.add(fileMenu);

		openSessionItem = new JMenuItem("Open Session", KeyEvent.VK_O);
		openSessionItem.setToolTipText("Open a saved Warp Drive session");
		openSessionItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				openSession_actionPerformed(e);
			}
		});

		saveSessionItem = new JMenuItem("Save Session", KeyEvent.VK_S);
		saveSessionItem.setToolTipText("Save this Warp Drive session");
		saveSessionItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveSession_actionPerformed(e);
			}
		});

		exitItem = new JMenuItem("Exit", KeyEvent.VK_X);
		exitItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				exit_actionPerformed(e);
			}
		});

		// TODO: delete these lines once the functionality for these buttons is
		// implemented
		openSessionItem.setEnabled(false);
		saveSessionItem.setEnabled(false);

		fileMenu.add(openSessionItem);
		fileMenu.add(saveSessionItem);
		fileMenu.addSeparator();
		fileMenu.add(exitItem);
		// end building file menu

		// build the tools menu
		toolsMenu = new JMenu("Tools");
		toolsMenu.setMnemonic(KeyEvent.VK_T);
		menuBar.add(toolsMenu);

		optionsItem = new JMenuItem("Options...", KeyEvent.VK_O);
		optionsItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				options_actionPerformed(e);
			}
		});

		toolsMenu.add(optionsItem);
		// end building edit menu

		this.setJMenuBar(menuBar);
	}

	public void addPainterContainer(PainterContainer pc) {
		pcs.add(pc);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == ok) {
			WarpOptions opts = pane.parseOptions();
			PainterContainer pc = null;
			for (int i = 0; i < pcs.size(); i++) {
				if (pcs.get(i).getGenome().getVersion().equals(opts.genome)) {
					pc = pcs.get(i);
					break;
				}
			}
			if (genomeString!=null && !pcs.isEmpty())
				pc = pcs.get(0);
			if (pc == null) {
				/*
				 * this is a bit of a hack to let us create an initial
				 * RegionFrame. In theory, we should 1) be able to handle
				 * different types of Frames 2) know whether it's even
				 * appropriate to create a RegionFrame
				 */
				new RegionFrame(opts);
			} else {
				opts = pane.parseAndDiff();
				if (genomeString!=null)
					opts.genomeString = genomeString;
				pc.addPaintersFromOpts(opts);
				// if (pc instanceof RegionPanel) {
				// RegionPanel rp = (RegionPanel)pc;
				// System.err.println("Setting Region\n");
				// rp.setRegion(rp.getRegion());
				// }
				// if (pc instanceof Component) {
				// System.err.println("Trying to force a repaint of " + pc);
				// ((Component)pc).repaint();
				// }
			}
			if (genomeString == null)
				pane.close();
			else
				pane.setClosed();
			this.dispose();
		} else if (e.getSource() == cancel) {
			pane.close();
			this.dispose();
		}
	}

	/**
	 * 
	 * @param e
	 */
	void openSession_actionPerformed(ActionEvent e) {
		// TODO Implement this method
	}

	/**
	 * 
	 * @param e
	 */
	void saveSession_actionPerformed(ActionEvent e) {
		// TODO Implement this method
	}

	/**
	 * Exit the program
	 * 
	 * @param e
	 */
	void exit_actionPerformed(ActionEvent e) {
		// TODO Check if any data is being viewed, if so prompt for confirmation
		// but if only the WarpOptionFrame is open then just exit
		boolean dataWindowsOpen = true;

		if (dataWindowsOpen) {
			int confirmResult = JOptionPane.showConfirmDialog(this, "Are you sure you want to exit Warp Drive?",
					"Confirm Exit", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);

			if (confirmResult == JOptionPane.NO_OPTION) {
				return;
			}
		}
		pane.close();
		try {
			Thread.sleep(400);
		} catch (Exception ex) {

		}
		System.exit(0);
	}

	/**
	 * Open the dialog to set preferences
	 * 
	 * @param e
	 */
	void options_actionPerformed(ActionEvent e) {
		new WarpOptionsDialog(this, this, pane.parseOptions());
	}

	/**
	 * Configure log4j
	 */
	public static void configureLogging() {
		ClassLoader loader = WarpOptionsFrame.class.getClassLoader();
		PropertyConfigurator.configure(loader.getResource("edu/mit/csail/cgs/utils/config/log4j.properties"));
	}

	/**
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		System.out.println("Starting WarpDrive (version: 2018-01-21)\n");
		try {
			WarpOptionsFrame.configureLogging();

			WarpOptions opts = WarpOptions.parseCL(args);
			new WarpOptionsFrame(opts);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
