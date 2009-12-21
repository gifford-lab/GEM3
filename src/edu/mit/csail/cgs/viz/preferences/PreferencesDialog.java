package edu.mit.csail.cgs.viz.preferences;

import java.util.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.viz.utils.GenomeSelectPanel;

public class PreferencesDialog extends JDialog {
	
	private PreferencesModel.Default model;
	private Map<String,JComponent> components;
	private PreferencesComponentTranslator translator;

	public PreferencesDialog(PreferencesModel.Default mod) {
		model = mod;
		translator = new DefaultComponentTranslator();
		components = new HashMap<String,JComponent>();
		config();
	}

	public PreferencesDialog(Frame p, PreferencesModel.Default mod) {
		super(p);
		model = mod;
		translator = new DefaultComponentTranslator();
		components = new HashMap<String,JComponent>();
		config();
	}

	private void config() { 
		Runnable r = new Runnable() { 
			public void run() { 

				Set<String> keys = model.getKeys();
				int numKeys = keys.size();

				Container cnt = (Container)getContentPane();
				cnt.setLayout(new BorderLayout());
				
				JPanel c = new JPanel();
				c.setLayout(new BorderLayout());
				
				JPanel center = new JPanel();
				//center.setLayout(new GridLayout(numKeys, 1));
				center.setLayout(new BoxLayout(center, BoxLayout.Y_AXIS));
				center.setBorder(new TitledBorder("Preferences"));
				
				for(String key : keys) { 
					Object val = model.getValue(key);
					JComponent comp = translator.getComponent(key, val);
					components.put(key, comp);
					center.add(comp);
				}
				
				c.add(center, BorderLayout.CENTER);
				
				JPanel buttons = new JPanel();
				buttons.setLayout(new FlowLayout());

				JButton ok = new JButton("Ok");
				JButton cancel = new JButton("Cancel");
				
				ok.addActionListener(new ActionListener() { 
					public void actionPerformed(ActionEvent e) { 
						doOK();
					}
				});
				cancel.addActionListener(new ActionListener() { 
					public void actionPerformed(ActionEvent e) { 
						doCancel();
					}
				});
				
				buttons.add(ok);
				buttons.add(cancel);
				
				c.add(buttons, BorderLayout.SOUTH);
				
				PreferencesDialog.this.addWindowListener(new WindowAdapter() {

					public void windowClosing(WindowEvent arg0) {
						System.out.println("Closing.");
						doCancel();
					}

				});
				
				cnt.add(new JScrollPane(c), BorderLayout.CENTER);
				
				setVisible(true);
				pack();
			}
		};
		
		EventQueue.invokeLater(r);
	}
	
	private void doOK() { 
		Runnable r = new Runnable() { 
			public void run() { 
				System.out.println("Performing OK.");
				setVisible(false);
				
				for(String k : components.keySet()) { 
					JComponent comp = components.get(k);
					Object val = translator.getValue(model.getValue(k), comp);
					model.setValue(k, val);
				}
				
				model.fireUpdatedEvent(PreferencesDialog.this);
				dispose();
			}
		};
		
		EventQueue.invokeLater(r);
	}
	
	private void doCancel() { 
		Runnable r = new Runnable() { 
			public void run() { 
				System.out.println("Performing CANCEL.");
				setVisible(false);
				model.fireCanceledEvent(PreferencesDialog.this);
				dispose();
			}
		};		
		
		EventQueue.invokeLater(r);
	}
}

interface PreferencesComponentTranslator { 
	public JComponent getComponent(String key, Object val);
	public Object getValue(Object oldValue, JComponent comp);
}

class LabeledComponent extends JPanel {
	
	private String label;
	private JComponent body;
	
	public LabeledComponent(String v, JComponent c) { 
		super();
		setLayout(new GridLayout(1, 2));
		label = v;
		body = c;
		JLabel lbl = new JLabel(label);
		lbl.setHorizontalAlignment(JLabel.RIGHT);
		add(lbl);
		add(body);
	}
	
	public JComponent getBody() { return body; }
}

class DefaultComponentTranslator implements PreferencesComponentTranslator {
	
	public DefaultComponentTranslator() {}
	
	public Object getValue(Object val, JComponent comp) { 
		if(val instanceof String) { 
			LabeledComponent lc = (LabeledComponent)comp;
			JTextField tf = (JTextField)lc.getBody();
			//JTextField tf = (JTextField)comp;
			return tf.getText();
		} 
		
		if(val instanceof Integer) { 
			LabeledComponent lc = (LabeledComponent)comp;
			JTextField tf = (JTextField)lc.getBody();
			//JTextField tf = (JTextField)comp;
			return Integer.parseInt(tf.getText());
		}
		
		if(val instanceof Genome) { 
			GenomeSelectPanel gsp = (GenomeSelectPanel)comp;
			try {
				return Organism.findGenome(gsp.getGenome());
			} catch (NotFoundException e) {
				e.printStackTrace();
				return val;
			}
		}
		
		if(val instanceof Boolean) { 
			JCheckBox cb = (JCheckBox)comp;
			return cb.isSelected();
		}
		
		if(val instanceof File) { 
			FileSelectPanel fsp = (FileSelectPanel)comp;
			return fsp.getFile();
		}

		return val;
	}
	
	public JComponent getComponent(String key, Object val) { 

		if(val instanceof String) { 
			return new LabeledComponent(key, new JTextField(val.toString()));
			//JComponent c = new JTextField(val.toString());
			//c.setBorder(new TitledBorder(key));
			//return c;
		} 
		
		if(val instanceof Integer) { 
			return new LabeledComponent(key, new JTextField(val.toString()));
			//JComponent c = new JTextField(val.toString());
			//c.setBorder(new TitledBorder(key));
			//return c;
		}
		
		if(val instanceof Genome) { 
			GenomeSelectPanel gsp = new GenomeSelectPanel((Genome)val);
			gsp.setBorder(new TitledBorder(key));
			return gsp;
		}
		
		if(val instanceof Boolean) { 
			JCheckBox cb = new JCheckBox(key, (Boolean)val);
			cb.setHorizontalTextPosition(JCheckBox.LEFT);
			return cb;
		}
		
		if(val instanceof File) { 
			FileSelectPanel fsp = new FileSelectPanel((File)val);
			fsp.setBorder(new TitledBorder(key));
			return fsp;
		}
		
		return new LabeledComponent(key, new JTextField(val.toString()));
	}
}
