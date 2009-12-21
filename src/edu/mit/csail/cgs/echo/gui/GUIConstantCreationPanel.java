/*
 * Created on Apr 12, 2007
 */
package edu.mit.csail.cgs.echo.gui;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.TitledBorder;

import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.echo.*;
import edu.mit.csail.cgs.ewok.types.SelfDescribingConstant;
import edu.mit.csail.cgs.ewok.types.ValueWrapper;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.viz.components.SelectionEvent;

public class GUIConstantCreationPanel extends JPanel implements EventSource<CreationEvent<SelfDescribingConstant>> {
    
    private DefaultComboBoxModel typeModel;
    private JComboBox typeBox;
    private JTextArea entryField;
    private JButton createButton;
    
    private EventSource.Default<CreationEvent<SelfDescribingConstant>> src;

    public GUIConstantCreationPanel() { 
        super();
        
        src = new EventSource.Default<CreationEvent<SelfDescribingConstant>>();
        
        typeModel = new DefaultComboBoxModel();
        typeBox = new JComboBox(typeModel);
        entryField = new JTextArea();
        createButton = new JButton("Create Constant");

        setLayout(new BorderLayout());
        
        JPanel typePanel = new JPanel(); 
        typePanel.setLayout(new BorderLayout());
        typePanel.setBorder(new TitledBorder("Constant Type:"));
        typePanel.add(typeBox, BorderLayout.NORTH);
        add(typePanel, BorderLayout.NORTH);
        
        JPanel entryPanel = new JPanel();
        entryPanel.setLayout(new BorderLayout());
        entryPanel.setBorder(new TitledBorder("Constant Type:"));
        entryPanel.add(entryField, BorderLayout.CENTER);
        add(entryPanel, BorderLayout.CENTER);

        add(createButton, BorderLayout.SOUTH);
        
        createButton.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                createConstant();
            }
        });
        

        registerStringFactory(new ConstantParser.StringParser());
        registerStringFactory(new ConstantParser.DoubleParser());
        registerStringFactory(new ConstantParser.IntegerParser());
        registerStringFactory(new ConstantParser.GenomeParser());
        registerFactory(new ExptLocatorFactory());
        registerFactory(new BindingScanFactory());
        registerStringFactory(new ConstantParser.AnnotationParamsParser());
        registerFactory(new FileFactory());
        registerStringFactory(new ConstantParser.RegionParser());
        
        typeModel.setSelectedItem(typeModel.getElementAt(0));
    }
    
    public void addEventListener(Listener<CreationEvent<SelfDescribingConstant>> el) {
        src.addEventListener(el);
    }

    public boolean hasListeners() {
        return src.hasListeners();
    }

    public void removeEventListener(Listener<CreationEvent<SelfDescribingConstant>> el) {
        src.removeEventListener(el);
    }
    
    public void registerFactory(SelfDescribingConstantFactory f) {
        System.out.println("Registering: " + f.getClass().getName());
        typeModel.addElement(f);
    }

    public void registerStringFactory(ConstantParser parser) {
        SelfDescribingConstantFactory fct = new EntryFieldWrapper(parser.toString(), parser);
        registerFactory(fct);
    }
    
    public void createConstant() { 
        Object sel = typeModel.getSelectedItem();
        System.out.println("Selected: " + sel.getClass().getName());
        SelfDescribingConstantFactory constFactory = (SelfDescribingConstantFactory)typeModel.getSelectedItem();
        if(constFactory != null) { 
            if(constFactory.isImmediate()) { 
                SelfDescribingConstant k = constFactory.createObject();
                src.fireEvent(new CreationEvent<SelfDescribingConstant>(this, k));
            } else { 
                DelayedCreationFactory dcf = new DelayedCreationFactory(src.getListeners(), constFactory);
                Thread t = new Thread(dcf);
                t.start();
            }
        }
    }
    
    private static class FileFactory implements SelfDescribingConstantFactory {
        
        public FileFactory() {}

        public boolean isImmediate() {
            return false;
        }

        public SelfDescribingConstant createObject() {
            JFileChooser chooser;
            chooser = new JFileChooser(new File(System.getProperty("user.dir")));
            int v = chooser.showSaveDialog(null);
            if(v == JFileChooser.APPROVE_OPTION) { 
                File f = chooser.getSelectedFile();
                return new ValueWrapper(f);
            } else { 
                return null;
            }
        }
        
        public String toString() { return "File"; }
    }
    
    private class BindingScanFactory 
    	implements SelfDescribingConstantFactory, Listener<CreationEvent<BindingScan>> {
    	
    	private BindingScan selected;
    	
    	public BindingScanFactory() { 
    		selected = null;
    	}

		public boolean isImmediate() {
			return false;
		}
		
		public String toString() { return "BindingScan"; }

		public SelfDescribingConstant createObject() {
            String text = entryField.getText().trim();
            Genome genome = null;
            
            try {
                genome = Organism.findGenome(text);
    			BindingSelectFrame f = new BindingSelectFrame(genome);
    			f.addEventListener(this);

            } catch (NotFoundException e1) {
                genome = null;
            } catch (SQLException e) {
				e.printStackTrace();
				return null;
			}
            
			synchronized(this) { 
				try {
					wait();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			return new ValueWrapper(selected);
		}

		public void eventRegistered(CreationEvent<BindingScan> e) {
            selected = e.getValue();
            synchronized(this) { 
                notifyAll();
            }
		}
    	
    }
    
    private class ExptLocatorFactory 
    	implements SelfDescribingConstantFactory, Listener<CreationEvent<ExptLocator>> {
    	
    	private ExptLocator selected;
    	
    	public ExptLocatorFactory() { 
    		selected = null;
    	}
    	
		public boolean isImmediate() {
			return false;
		}
		
		public String toString() { return "ExptLocator"; }

		public SelfDescribingConstant createObject() {
            String text = entryField.getText().trim();
            Genome genome = null;
            
            try {
                genome = Organism.findGenome(text);
            } catch (NotFoundException e1) {
                genome = null;
            }
            
			ExptSelectFrame f = new ExptSelectFrame(genome);
			f.addEventListener(this);
			synchronized(this) { 
				try {
					wait();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			return new ValueWrapper(selected);
		}

		public void eventRegistered(CreationEvent<ExptLocator> e) {
            selected = e.getValue();
            synchronized(this) { 
                notifyAll();
            }
		} 
    }
    
    private class EntryFieldWrapper implements SelfDescribingConstantFactory { 

        private String name;
        private ConstantParser parser;
        
        public EntryFieldWrapper(String n, ConstantParser p) {
            name = n;
            parser = p;
        }
        
        public String toString() { return name; }
        
        public boolean isImmediate() { return true; }

        public SelfDescribingConstant createObject() {
            String text = entryField.getText().trim();
            SelfDescribingConstant k = parser.parseConstant(text);
            entryField.setText("");
            return k;
        }
    }

}
