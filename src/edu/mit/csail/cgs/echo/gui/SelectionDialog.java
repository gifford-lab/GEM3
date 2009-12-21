package edu.mit.csail.cgs.echo.gui;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.types.ValueWrapper;
import edu.mit.csail.cgs.utils.NotFoundException;

import java.util.*;

public class SelectionDialog<X> extends JDialog {
	
	public static final int OK = 0;
	public static final int CANCEL = 1;
	
	private SwingSelectionComponent<X> selComp;
	private JButton ok, cancel;
	private int choice;
	private EchoAction<X> okAction;
    
    private LinkedList<SelectionDialogListener> listeners;

	public SelectionDialog(SwingSelectionComponent<X> comp, EchoAction<X> okact) { 
		super();

		selComp = comp;
		okAction = okact;
        
        listeners = new LinkedList<SelectionDialogListener>();
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		
		JPanel buttonPanel = new JPanel(); buttonPanel.setLayout(new FlowLayout());
		buttonPanel.add(ok = new JButton("Ok"));
		buttonPanel.add(cancel = new JButton("Cancel"));
		
		c.add(buttonPanel, BorderLayout.SOUTH);
		c.add(comp.asJComponent(), BorderLayout.CENTER);
		
		choice = CANCEL;
		
		ok.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				choice = OK;
				okAction.doAction(selComp.getSelectedValue());
                closeDialog();
			}
		});
        
		cancel.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				choice = CANCEL;
                closeDialog();
			}
		});		
	}

    
    public void addSelectionDialogListener(SelectionDialogListener l) { listeners.addLast(l); }
    public void removeSelectionDialogListener(SelectionDialogListener l) { listeners.remove(l); }
    
    private void fireSelectionDialogClosed() { 
        SelectionDialogEvent evt = new SelectionDialogEvent(this);
        for(SelectionDialogListener sl : listeners) { 
            sl.selectionDialogClosed(evt);
        }
    }
    
    public void openDialog() { 
        setVisible(true);
        pack();
    }
    
    public void closeDialog() { 
        setVisible(false);
        fireSelectionDialogClosed();
    }

	public X getSelectedValue() { return selComp.getSelectedValue(); }
	public boolean isAcceptedChoice() { return choice == OK; }
    
    private static class WrapperAction<X> implements EchoAction<X> { 
        protected ValueWrapper wrapper;
        public WrapperAction(ValueWrapper vw) { wrapper = vw; }
        
        public void doAction(X param) {
            wrapper.setConstantValue(param);
        }
    }

    public static SelectionDialog fromValueWrapper(ValueWrapper wrapper) {
                
        if(wrapper.getConstantClass().equals(String.class)) { 
            SwingSelectionComponent<String> ssc = new StringSelectionPanel((String)wrapper.getConstantValue());
            EchoAction<String> act = new WrapperAction<String>(wrapper);
            return new SelectionDialog<String>(ssc, act);
        }

        if(wrapper.getConstantClass().equals(Double.class)) { 
            SwingSelectionComponent<String> ssc = new StringSelectionPanel(String.valueOf((Double)wrapper.getConstantValue()));
            EchoAction<String> act = new WrapperAction<String>(wrapper);
            return new SelectionDialog<String>(ssc, act);
        }

        if(wrapper.getConstantClass().equals(Integer.class)) { 
            SwingSelectionComponent<String> ssc = new StringSelectionPanel(String.valueOf((Integer)wrapper.getConstantValue()));
            EchoAction<String> act = new WrapperAction<String>(wrapper);
            return new SelectionDialog<String>(ssc, act);
        }
        
        if(wrapper.getConstantClass().equals(Genome.class)) { 
            SwingSelectionComponent<String> ssc = new StringSelectionPanel(((Genome)wrapper.getConstantValue()).getVersion());
            EchoAction<String> act = new WrapperAction<String>(wrapper) { 
                public void doAction(String val) { 
                    try {
                        Genome newGenome = Organism.findGenome(val);
                        wrapper.setConstantValue(newGenome);
                    } catch (NotFoundException e) {
                        e.printStackTrace();
                    }
                }
            };
            return new SelectionDialog<String>(ssc, act);            
        }

        return null;
    }


}
