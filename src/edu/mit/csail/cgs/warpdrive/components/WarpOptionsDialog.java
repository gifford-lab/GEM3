package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.*;
import java.beans.*;

import edu.mit.csail.cgs.warpdrive.*;
import edu.mit.csail.cgs.utils.GenericFileFilter;

import org.apache.log4j.Logger;

/**
 * This is a dialog for managing the options/preferences for the WarpOptions
 * GUI. Unfortunately since the main GUI for Warp Drive is called 
 * WarpOptionsFrame the name is a bit confusing.
 *  
 * @author Bob
 *
 */
public class WarpOptionsDialog extends JDialog {

	static Logger logger = Logger.getLogger(WarpOptionsDialog.class);
	
	private WarpOptionsFrame optionsFrame;
	private WarpOptions options;
	
	private JPanel mainPanel = new JPanel();
	private GridBagLayout mainPanelLayout = new GridBagLayout();
	
	//gui components for various preferences
	private JTabbedPane tabbedPane = new JTabbedPane();
	
	private JPanel generalTab = new JPanel();
	private GridBagLayout generalTabLayout = new GridBagLayout();
	
	private JPanel importExportPanel = new JPanel();
	private GridBagLayout importExportPanelLayout = new GridBagLayout();
	private Border importExportPanelBorder;
	private JButton exportOptionsButton = new JButton("Save and Export");
	private JButton importOptionsButton = new JButton("Import and Load");
	private JButton originalsButton = new JButton("Restore Original Settings");

	private JPanel windowAppearancePanel = new JPanel();
	private GridBagLayout windowAppearancePanelLayout = new GridBagLayout();
	private Border windowAppearancePanelBorder;
	private JPanel windowSizePanel = new JPanel();
	private GridBagLayout windowSizePanelLayout = new GridBagLayout();
	private JLabel mainWindowSizeLabel = new JLabel("Main Window Size:");
	private DecimalFormat mainWindowWidthFieldFormat = new DecimalFormat("0");
	private JFormattedTextField mainWindowWidthField = new JFormattedTextField(mainWindowWidthFieldFormat);
	private JLabel mainWindowSizeXLabel = new JLabel("x");
	private DecimalFormat mainWindowHeightFieldFormat = new DecimalFormat("0");
	private JFormattedTextField mainWindowHeightField = new JFormattedTextField(mainWindowHeightFieldFormat);
	
	private JPanel windowLocationPanel = new JPanel();
	private GridBagLayout windowLocationPanelLayout = new GridBagLayout();
	private JLabel mainWindowLocationLabel = new JLabel("Main Window Location:");
	private ButtonGroup locationButtonGroup = new ButtonGroup();
	private JRadioButton centeredButton = new JRadioButton("Centered on Screen");
	private JRadioButton topLeftButton = new JRadioButton("Top Left Corner at X:");
	private DecimalFormat topLeftXFieldFormat = new DecimalFormat("0");
	private JFormattedTextField topLeftXField = new JFormattedTextField(topLeftXFieldFormat);
	private JLabel topLeftYLabel = new JLabel("Y:");
	private DecimalFormat topLeftYFieldFormat = new DecimalFormat("0");
	private JFormattedTextField topLeftYField = new JFormattedTextField(topLeftYFieldFormat); 
	
	private JPanel chipSeqTab = new JPanel();
	
	
	//panel for ok/cancel buttons
	private JPanel buttonPanel = new JPanel();
	private GridBagLayout buttonPanelLayout = new GridBagLayout();
	private JCheckBox saveAsDefaultsCB = new JCheckBox("Save As Defaults");
	private JButton saveButton = new JButton("Save Options");
	private JButton cancelButton = new JButton("Cancel");
	
	/**
	 * Construct the dialog box
	 * @param owner This will typically be the WarpOptionsFrame
	 */
	public WarpOptionsDialog(JFrame owner, WarpOptionsFrame wof, WarpOptions options) {
		super(owner, true);
		this.optionsFrame = wof;
		this.options = options;
		
		try {
			init();
		}
		catch(Exception ex) {
			logger.error("Exception while initializing WarpOptionsDialog", ex);
		}
	}
	
	
	/**
	 * Initialize the GUI
	 * @throws Exception if an exception occurs during initialization
	 */
	private void init() throws Exception {
		this.setTitle("Warp Drive Options");
	    this.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
	    this.getContentPane().add(mainPanel, BorderLayout.CENTER);

	    mainPanel.setLayout(mainPanelLayout);

	    //initialize tabbed pane
	    
	    //initialize general tab
	    generalTab.setLayout(generalTabLayout);
	    importExportPanel.setLayout(importExportPanelLayout);
	    importExportPanelBorder = BorderFactory.createCompoundBorder(new TitledBorder(BorderFactory.createEtchedBorder(Color.
	            white, new Color(165, 163, 151)), "Warp Drive Options:"), BorderFactory.createEmptyBorder(6, 6, 5, 5));
	    importExportPanel.setBorder(importExportPanelBorder);

	    exportOptionsButton.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	          exportOptionsButton_actionPerformed(e);
	        }
	      });
	    
	    importOptionsButton.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	          importOptionsButton_actionPerformed(e);
	        }
	      });
	    
	    originalsButton.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	          originalsButton_actionPerformed(e);
	        }
	      });
	    
	    
	    importExportPanel.add(exportOptionsButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
	    importExportPanel.add(importOptionsButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
	    importExportPanel.add(originalsButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
	    importExportPanel.add(Box.createVerticalStrut(0), new GridBagConstraints(3, 0, 1, GridBagConstraints.REMAINDER, 1.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
	    
	    windowAppearancePanel.setLayout(windowAppearancePanelLayout);
	    windowAppearancePanelBorder = BorderFactory.createCompoundBorder(new TitledBorder(BorderFactory.createEtchedBorder(Color.
	            white, new Color(165, 163, 151)), "Main Window Appearance:"), BorderFactory.createEmptyBorder(6, 6, 5, 5));
	    windowAppearancePanel.setBorder(windowAppearancePanelBorder);
	    
	    windowSizePanel.setLayout(windowSizePanelLayout);
	    
	    mainWindowWidthField.setColumns(5);
	    mainWindowWidthField.setHorizontalAlignment(JTextField.RIGHT);
	    mainWindowWidthField.setText(mainWindowWidthFieldFormat.format(optionsFrame.getWidth()));
	    mainWindowWidthField.addPropertyChangeListener("value", new PropertyChangeListener() {
	    	public void propertyChange(PropertyChangeEvent e) {
	    		mainWindowWidth_propertyChange(e);
	    	}
	    });

	    mainWindowHeightField.setColumns(5);
	    mainWindowHeightField.setHorizontalAlignment(JTextField.RIGHT);
	    mainWindowHeightField.setText(mainWindowHeightFieldFormat.format(optionsFrame.getHeight()));
	    mainWindowHeightField.addPropertyChangeListener("value", new PropertyChangeListener() {
	    	public void propertyChange(PropertyChangeEvent e) {
	    		mainWindowHeight_propertyChange(e);
	    	}
	    });

	    windowSizePanel.add(mainWindowSizeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
	    windowSizePanel.add(mainWindowWidthField, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
	    windowSizePanel.add(mainWindowSizeXLabel, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 0), 0, 0));
	    windowSizePanel.add(mainWindowHeightField, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 2, 0, 0), 0, 0));
	    windowSizePanel.add(Box.createVerticalStrut(0), new GridBagConstraints(4, 0, 1, GridBagConstraints.REMAINDER, 1.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
	    
	    windowLocationPanel.setLayout(windowLocationPanelLayout);
	    
	    locationButtonGroup.add(centeredButton);
	    locationButtonGroup.add(topLeftButton);
	    
	    topLeftButton.addItemListener(new ItemListener() {
	    	public void itemStateChanged(ItemEvent e) {
	    		topLeftButton_stateChanged(e);
	    	}
	    });
	    
	    topLeftXField.setColumns(5);
	    topLeftXField.setHorizontalAlignment(JTextField.RIGHT);
	    topLeftXField.setText(topLeftXFieldFormat.format(optionsFrame.getX()));
	    topLeftXField.addPropertyChangeListener("value", new PropertyChangeListener() {
	    	public void propertyChange(PropertyChangeEvent e) {
	    		topLeftX_propertyChange(e);
	    	}
	    });
	    
	    topLeftYField.setColumns(5);
	    topLeftYField.setHorizontalAlignment(JTextField.RIGHT);
	    topLeftYField.setText(topLeftYFieldFormat.format(optionsFrame.getY()));
	    topLeftYField.addPropertyChangeListener("value", new PropertyChangeListener() {
	    	public void propertyChange(PropertyChangeEvent e) {
	    		topLeftY_propertyChange(e);
	    	}
	    });

	    
	    if (options.isWindowCentered()) {
	    	centeredButton.setSelected(true);
	    	topLeftXField.setEnabled(false);
	    	topLeftYField.setEnabled(false);
	    }
	    else {
	    	topLeftButton.setSelected(true);
	    }
	    
	    windowLocationPanel.add(mainWindowLocationLabel, new GridBagConstraints(0, 0, 4, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
	    windowLocationPanel.add(centeredButton, new GridBagConstraints(0, 1, 4, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 18, 0, 0), 0, 0));
	    windowLocationPanel.add(topLeftButton, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 18, 0, 0), 0, 0));
	    windowLocationPanel.add(topLeftXField, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
	    windowLocationPanel.add(topLeftYLabel, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
	    windowLocationPanel.add(topLeftYField, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 6, 0, 0), 0, 0));
	    windowLocationPanel.add(Box.createVerticalStrut(0), new GridBagConstraints(4, 0, 1, GridBagConstraints.REMAINDER, 1.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
	    
	    
	    windowAppearancePanel.add(windowSizePanel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 6, 0, 5), 0, 0));
	    windowAppearancePanel.add(windowLocationPanel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(6, 6, 5, 5), 0, 0));
	    windowAppearancePanel.add(Box.createVerticalStrut(0), new GridBagConstraints(1, 0, 1, GridBagConstraints.REMAINDER, 1.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
	    
	    generalTab.add(importExportPanel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 6, 0, 5), 0, 0));
	    generalTab.add(windowAppearancePanel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
	            , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(6, 6, 5, 5), 0, 0));
	    generalTab.add(Box.createHorizontalStrut(0), new GridBagConstraints(0, 2, GridBagConstraints.REMAINDER, 1, 0.0, 1.0
	            , GridBagConstraints.WEST, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
	    //done initializing general options tab
	    
	    
	    tabbedPane.add("General", generalTab);
	    tabbedPane.add("Chip-Seq", chipSeqTab);
	    //done initializing tabbedPane

	    
	    //initialize button panel
	    saveAsDefaultsCB.setSelected(true);
	    saveAsDefaultsCB.setToolTipText("Save these settings as defaults for future uses of Warp Drive.");
	    
	    buttonPanel.setLayout(buttonPanelLayout);
	    buttonPanel.add(Box.createVerticalStrut(0), new GridBagConstraints(0, 0, 1, GridBagConstraints.REMAINDER, 1.0, 0.0
	            , GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
	    buttonPanel.add(saveAsDefaultsCB, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
	    buttonPanel.add(saveButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
	    buttonPanel.add(cancelButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 6, 0, 0), 0, 0));
	        
	    saveButton.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	          saveButton_actionPerformed(e);
	        }
	      });

	    cancelButton.addActionListener(new ActionListener() {
	        public void actionPerformed(ActionEvent e) {
	          cancelButton_actionPerformed(e);
	        }
	      });

	    mainPanel.add(tabbedPane, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
	            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(12, 12, 0, 11), 0, 0));

	    mainPanel.add(buttonPanel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
	            , GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(12, 12, 11, 11), 0, 0));
	    
	    this.pack();
	    this.setLocationRelativeTo(getParent());
	    this.setResizable(false);
	    this.setVisible(true);
	}

	
	/**
	 * Check that all the option values are valid
	 * @return true if all the option values are valid, false if any are not
	 */
	private boolean checkOptionValues() {
		int width = Integer.valueOf(mainWindowWidthField.getText()).intValue();
		if (!options.checkPreferredWindowWidth(width)) {
			String errorMessage = "Main Window width must be between " + WarpOptions.MIN_WINDOW_WIDTH
			+ " and " + WarpOptions.MAX_WINDOW_WIDTH + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Width",
                    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		int height = Integer.valueOf(mainWindowHeightField.getText()).intValue();
		if (!options.checkPreferredWindowHeight(height)) {
			String errorMessage = "Main Window height must be between " + WarpOptions.MIN_WINDOW_HEIGHT
			+ " and " + WarpOptions.MAX_WINDOW_HEIGHT + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Height",
                    JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		boolean isCentered = centeredButton.isSelected();
		if (!isCentered) {
			int topLeftX = Integer.valueOf(topLeftXField.getText()).intValue();
			if (!options.checkPreferredTopLeftX(topLeftX)) {
				String errorMessage = "Top Left X location must be between " + WarpOptions.MIN_TOP_LEFT_X
				+ " and " + WarpOptions.MAX_TOP_LEFT_X + ".";
				JOptionPane.showMessageDialog(this,
	                    errorMessage,
	                    "Invalid Top Left X Location",
	                    JOptionPane.ERROR_MESSAGE);
				return false;
			}
			
			int topLeftY = Integer.valueOf(topLeftYField.getText()).intValue();
			if (!options.checkPreferredTopLeftY(topLeftY)) {
				String errorMessage = "Top Left Y location must be between " + WarpOptions.MIN_TOP_LEFT_Y
				+ " and " + WarpOptions.MAX_TOP_LEFT_Y + ".";
				JOptionPane.showMessageDialog(this,
	                    errorMessage,
	                    "Invalid Top Left Y Location",
	                    JOptionPane.ERROR_MESSAGE);
				return false;
			}
		}
		
		//if all checks have passed then return true;
		return true;
	}
	
		
	/**
	 * Save the options to the WarpOptionsModel and optionally set them as 
	 * defaults for Warp Drive
	 * @param saveAsDefaults
	 * @return whether the save succeeded
	 */	
	private boolean saveOptions(boolean saveAsDefaults) {
		if (checkOptionValues()) {
			//get the settings from the text fields
			int newWidth = Integer.valueOf(mainWindowWidthField.getText());
			int newHeight = Integer.valueOf(mainWindowHeightField.getText());
			boolean isCentered = centeredButton.isSelected();
			int newTopLeftX = Integer.valueOf(topLeftXField.getText());
			int newTopLeftY = Integer.valueOf(topLeftYField.getText());

			//all of these sets should work because checkOptionValues worked
			boolean valid = true;
			valid = valid && options.setPreferredWindowWidth(newWidth);
			valid = valid && options.setPreferredWindowHeight(newHeight);
			options.setWindowCentered(isCentered);
			if (!isCentered) {
				valid = valid && options.setPreferredTopLeftX(newTopLeftX);
				valid = valid && options.setPreferredTopLeftY(newTopLeftY);
			}
			
			if (valid) {
				if (saveAsDefaults) {
					options.saveOptions();
				}
			}
			else {
				String errorMessage = "Error saving options. See log for details.";
				JOptionPane.showMessageDialog(this,
						errorMessage,
						"Error saving options",
						JOptionPane.ERROR_MESSAGE);
				logger.error("Error saving options, checkOptionValues passed, but a set method returned false.");
			}
			return valid;
		}
		else {
			return false;
		}
	}
		
		
	/**
	 * Respond to changes in the value in the mainWindowWidthField. Check that the
	 * new value is valid and if not then revert to a valid value.
	 * 
	 * @param e PropertyChangeEvent
	 */
	void mainWindowWidth_propertyChange(PropertyChangeEvent e) {
		double newValue;
		if (e.getNewValue() instanceof Long) {
			Long newValueLong = (Long)e.getNewValue();
			newValue = newValueLong.doubleValue();
		}
		else if (e.getNewValue()instanceof Double) {
			Double newValueDouble = (Double)e.getNewValue();
			newValue = newValueDouble.doubleValue();
		}
		else {
			//new value is some unusual type, so just let it go through
			return;
		}
		//figure out what the text value will be and check that it's valid
		int textValue = Long.valueOf(mainWindowWidthFieldFormat.format(newValue)).intValue();
		if ((newValue > WarpOptions.MAX_WINDOW_WIDTH) || (newValue < WarpOptions.MIN_WINDOW_WIDTH)
				|| !options.checkPreferredWindowWidth(textValue)) {
			String errorMessage = "Main Window width must be between " + WarpOptions.MIN_WINDOW_WIDTH
			+ " and " + WarpOptions.MAX_WINDOW_WIDTH + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Width",
                    JOptionPane.ERROR_MESSAGE);
			
			//reset the value to something valid
			if (e.getOldValue() != null) {
				mainWindowWidthField.setValue(e.getOldValue());
			}
			else {
				mainWindowWidthField.setValue(new Integer(WarpOptions.DEFAULT_WINDOW_WIDTH));
			}
		}
	}

	
	/**
	 * Respond to changes in the value in the mainWindowHeightField. Check that the
	 * new value is valid and if not then revert to a valid value.
	 * 
	 * @param e PropertyChangeEvent
	 */
	void mainWindowHeight_propertyChange(PropertyChangeEvent e) {
		double newValue;
		if (e.getNewValue() instanceof Long) {
			Long newValueLong = (Long)e.getNewValue();
			newValue = newValueLong.doubleValue();
		}
		else if (e.getNewValue()instanceof Double) {
			Double newValueDouble = (Double)e.getNewValue();
			newValue = newValueDouble.doubleValue();
		}
		else {
			//new value is some unusual type, so just let it go through
			return;
		}
		//figure out what the text value will be and check that it's valid
		int textValue = Long.valueOf(mainWindowHeightFieldFormat.format(newValue)).intValue();
		if ((newValue > WarpOptions.MAX_WINDOW_HEIGHT) || (newValue < WarpOptions.MIN_WINDOW_HEIGHT)
				|| !options.checkPreferredWindowHeight(textValue)) {
			String errorMessage = "Main Window height must be between " + WarpOptions.MIN_WINDOW_HEIGHT
			+ " and " + WarpOptions.MAX_WINDOW_HEIGHT + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Height",
                    JOptionPane.ERROR_MESSAGE);
			
			//reset the value to something valid
			if (e.getOldValue() != null) {
				mainWindowHeightField.setValue(e.getOldValue());
			}
			else {
				mainWindowHeightField.setValue(new Integer(WarpOptions.DEFAULT_WINDOW_HEIGHT));
			}
		}
	}

	
	/**
	 * Respond to changes in the value in the topLeftXField. Check that the
	 * new value is valid and if not then revert to a valid value.
	 * 
	 * @param e PropertyChangeEvent
	 */
	void topLeftX_propertyChange(PropertyChangeEvent e) {
		double newValue;
		if (e.getNewValue() instanceof Long) {
			Long newValueLong = (Long)e.getNewValue();
			newValue = newValueLong.doubleValue();
		}
		else if (e.getNewValue()instanceof Double) {
			Double newValueDouble = (Double)e.getNewValue();
			newValue = newValueDouble.doubleValue();
		}
		else {
			//new value is some unusual type, so just let it go through
			return;
		}
		//figure out what the text value will be and check that it's valid
		int textValue = Long.valueOf(topLeftXFieldFormat.format(newValue)).intValue();
		if ((newValue > WarpOptions.MAX_TOP_LEFT_X) || (newValue < WarpOptions.MIN_TOP_LEFT_X)
				|| !options.checkPreferredTopLeftX(textValue)) {
			String errorMessage = "Main Window top left X location must be between " 
				+ WarpOptions.MIN_TOP_LEFT_X + " and " + WarpOptions.MAX_TOP_LEFT_X + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Top Left X Location",
                    JOptionPane.ERROR_MESSAGE);
			
			//reset the value to something valid
			if (e.getOldValue() != null) {
				topLeftXField.setValue(e.getOldValue());
			}
			else {
				topLeftXField.setValue(new Integer(WarpOptions.DEFAULT_TOP_LEFT_X));
			}
		}
	}
	

	/**
	 * Respond to changes in the value in the topLeftYField. Check that the
	 * new value is valid and if not then revert to a valid value.
	 * 
	 * @param e PropertyChangeEvent
	 */
	void topLeftY_propertyChange(PropertyChangeEvent e) {
		double newValue;
		if (e.getNewValue() instanceof Long) {
			Long newValueLong = (Long)e.getNewValue();
			newValue = newValueLong.doubleValue();
		}
		else if (e.getNewValue()instanceof Double) {
			Double newValueDouble = (Double)e.getNewValue();
			newValue = newValueDouble.doubleValue();
		}
		else {
			//new value is some unusual type, so just let it go through
			return;
		}
		//figure out what the text value will be and check that it's valid
		int textValue = Long.valueOf(topLeftYFieldFormat.format(newValue)).intValue();
		if ((newValue > WarpOptions.MAX_TOP_LEFT_Y) || (newValue < WarpOptions.MIN_TOP_LEFT_Y)
				|| !options.checkPreferredTopLeftY(textValue)) {
			String errorMessage = "Main Window top left Y location must be between " 
				+ WarpOptions.MIN_TOP_LEFT_Y + " and " + WarpOptions.MAX_TOP_LEFT_Y + ".";
			JOptionPane.showMessageDialog(this,
                    errorMessage,
                    "Invalid Main Window Top Left Y Location",
                    JOptionPane.ERROR_MESSAGE);
			
			//reset the value to something valid
			if (e.getOldValue() != null) {
				topLeftYField.setValue(e.getOldValue());
			}
			else {
				topLeftYField.setValue(new Integer(WarpOptions.DEFAULT_TOP_LEFT_Y));
			}
		}
	}
	

	/**
	 * Save the current options to the underlying WarpOptions and Export them
	 * to a file
	 *
	 * @param e ActionEvent the event created by the user clicking on the button
	 */
	void exportOptionsButton_actionPerformed(ActionEvent e) {
		//Create a file chooser to let the user specify the file to export
		//Add a file filter for XML files
		JFileChooser exportFileChooser = new JFileChooser();
		String[] allowedExtensions = { "xml" };
		String description = "XML files";
		FileFilter filter = new GenericFileFilter(allowedExtensions, description);
		exportFileChooser.addChoosableFileFilter(filter);
		exportFileChooser.setMultiSelectionEnabled(false);
		int returnVal = exportFileChooser.showSaveDialog(this);
		//Only proceed if the user selected a file and hit OK
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			//clean up the filename if necessary, including enforcing a .xml extension
			String filename = exportFileChooser.getSelectedFile().getAbsolutePath();
			if (!filter.accept(exportFileChooser.getSelectedFile())) {
				if (filename.charAt(filename.length() - 1) == '.') {
					filename = filename + allowedExtensions[0];
				}
				else {
					filename = filename + "." + allowedExtensions[0];
				}
			}

			//Create a file object with the specified name, check if the file
			//already exists, and if so confirm that it's OK to overwrite
			File exportFile = null;
			try {
				exportFile = new File(filename);
			}
			catch (Exception ex) {
				JOptionPane.showMessageDialog(this, filename + " is not a valid filename.", 
						"Invalid Filename", JOptionPane.ERROR_MESSAGE);
				return;
			}

			if (exportFile.exists()) {
				int choice = JOptionPane.showConfirmDialog(this,
						filename + " already exists.\nDo you want to replace it?", "Confirm Save",
						JOptionPane.YES_NO_OPTION, JOptionPane.WARNING_MESSAGE);
				if (choice != JOptionPane.YES_OPTION) {
					return;
				}
			}

			//Open an output stream for the file object and execute the export
			FileOutputStream fos;
			try {
				//open a file output stream
				fos = new FileOutputStream(exportFile);

				//save the specified settings to the underlying WarpOptions and export
				//Note: this requires that the options be saved as defaults
				boolean saveSuccessful = this.saveOptions(true);  

				if (saveSuccessful) {
					try {
						options.exportOptions(fos);
					}
					catch (IOException ioex) {
						String errorMessage = "Error exporting options. See log for details.";
						JOptionPane.showMessageDialog(this,
								errorMessage,
								"Error saving options",
								JOptionPane.ERROR_MESSAGE);
						logger.error("Error saving options", ioex);
					}
					catch (WarpDriveException wdex) {
						String errorMessage = "Error saving options. See log for details.";
						JOptionPane.showMessageDialog(this,
								errorMessage,
								"Error saving options",
								JOptionPane.ERROR_MESSAGE);
						logger.error("Error exporting options", wdex);
					}
					finally {
						try {
							fos.close();
						} 
						catch (IOException ex) {
							logger.error("Error closing options file", ex);
						}
					}
				}
			}
			catch (FileNotFoundException fnfex) {
				String errorMessage = "Error saving options. See log for details.";
				JOptionPane.showMessageDialog(this,
						errorMessage,
						"Error saving options",
						JOptionPane.ERROR_MESSAGE);
				logger.error("Error exporting options", fnfex);
			}
		}
	}
	

	/**
	 * Import options from a file specified by the user
	 *
	 * @param e ActionEvent the event created by the user clicking on the button
	 */
	void importOptionsButton_actionPerformed(ActionEvent e) {
		//Create a file chooser to let the user specify the file to import
		//Add a file filter for XML files
		JFileChooser importFileChooser = new JFileChooser();
		String[] allowedExtensions = { "xml" };
		String description = "XML files";
		FileFilter filter = new GenericFileFilter(allowedExtensions, description);
		importFileChooser.addChoosableFileFilter(filter);
		importFileChooser.setMultiSelectionEnabled(false);
		int returnVal = importFileChooser.showOpenDialog(this);
		//Only proceed if the user selected a file and hit OK
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String filename = importFileChooser.getSelectedFile().getAbsolutePath();

			//Create a file object with the specified name, check if the file
			//already exists, and if so confirm that it's OK to overwrite
			File importFile = null;
			try {
				importFile = new File(filename);
			}
			catch (Exception ex) {
				JOptionPane.showMessageDialog(this, filename + " is not a valid filename.", 
						"Invalid Filename", JOptionPane.ERROR_MESSAGE);
				return;
			}

			//Open an input stream for the file object and execute the import
			FileInputStream fis;
			try {
				//open a file input stream
				fis = new FileInputStream(importFile);

				try {
					options.importOptions(fis);
				}
				catch (IOException ioex) {
					String errorMessage = "Error importing options. See log for details.";
					JOptionPane.showMessageDialog(this,
							errorMessage,
							"Error importing options",
							JOptionPane.ERROR_MESSAGE);
					logger.error("Error importing options", ioex);
				}
				catch (WarpDriveException wdex) {
					String errorMessage = "Error importing options. See log for details.";
					JOptionPane.showMessageDialog(this,
							errorMessage,
							"Error importing options",
							JOptionPane.ERROR_MESSAGE);
					logger.error("Error importing options", wdex);
				}
				finally {
					try {
						fis.close();
					} 
					catch (IOException ex) {
						logger.error("Error closing options file", ex);
					}
				}
			}
			catch (FileNotFoundException fnfex) {
				String errorMessage = "Error importing options. See log for details.";
				JOptionPane.showMessageDialog(this,
						errorMessage,
						"Error importing options",
						JOptionPane.ERROR_MESSAGE);
				logger.error("Error importing options", fnfex);
			}
		}
	}

	
	/**
	 * Restore the options to their original values - load the original options
	 * but don't save them as defaults until the user actually hits the save
	 * button
	 *
	 * @param e ActionEvent the event created by the user clicking on the button
	 */
	void originalsButton_actionPerformed(ActionEvent e) {
		mainWindowWidthField.setValue(new Integer(WarpOptions.DEFAULT_WINDOW_WIDTH));
		mainWindowHeightField.setValue(new Integer(WarpOptions.DEFAULT_WINDOW_HEIGHT));
		if (WarpOptions.DEFAULT_WINDOW_IS_CENTERED) {
			centeredButton.setSelected(true);
		}
		else {
			topLeftButton.setSelected(true);
			topLeftXField.setValue(new Integer(WarpOptions.DEFAULT_TOP_LEFT_X));
			topLeftYField.setValue(new Integer(WarpOptions.DEFAULT_TOP_LEFT_Y));
		}
	}
	
	
	/**
	 * Enable/Disable the fields for editing the location of the top left
	 * corner of the window in response to a change in the manner by which
	 * the window location is specified.
	 * 
	 * @param e ItemEvent the event created by the user clicking on the radio buttons
	 */
	void topLeftButton_stateChanged(ItemEvent e) {
		if (e.getStateChange() == ItemEvent.SELECTED) {
			topLeftXField.setEnabled(true);
			topLeftYField.setEnabled(true);
		}
		else {
			topLeftXField.setEnabled(false);
			topLeftYField.setEnabled(false);
		}
	}

	
	/**
	 * Save changes and close dialog
	 *
	 * @param e ActionEvent the event created by the user clicking on the button
	 */
	void saveButton_actionPerformed(ActionEvent e) {	
		boolean saveSuccessful = this.saveOptions(saveAsDefaultsCB.isSelected());
		if (saveSuccessful) {
			this.dispose();
		}
		else {
			String errorMessage = "Error saving options. See log for details.";
			JOptionPane.showMessageDialog(this,
					errorMessage,
					"Error saving options",
					JOptionPane.ERROR_MESSAGE);
		}
	}



	/**
	 * Abort changes and close dialog
	 *
	 * @param e ActionEvent the event created by the user clicking on the button
	 */
	void cancelButton_actionPerformed(ActionEvent e) {
		this.dispose();
	}
}
