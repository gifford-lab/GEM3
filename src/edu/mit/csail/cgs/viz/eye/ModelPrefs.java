/*
 * Author: tdanford
 * Date: Jan 26, 2009
 */
package edu.mit.csail.cgs.viz.eye;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.models.*;

public class ModelPrefs<T extends Model> extends JFrame {
	
	private PrefsPanel<T> panel;
	private JButton ok, cancel;
	private LinkedList<ModelListener<T>> listeners;
	
	public ModelPrefs(T mdl) { 
		super("Model");
		
		Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
		listeners = new LinkedList<ModelListener<T>>();
		
		panel = new PrefsPanel<T>(mdl);
		c.add(panel, BorderLayout.CENTER);
		panel.setPreferredSize(new Dimension(450, panel.numFields() * 50));
		
		JPanel buttons = new JPanel();
		buttons.setLayout(new FlowLayout());
		buttons.add(ok = new JButton(createOkAction()));
		buttons.add(cancel = new JButton(createCancelAction()));
		c.add(buttons, BorderLayout.SOUTH);
		
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	public void addModelListener(ModelListener<T> lst) { 
		listeners.add(lst);
	}
	
	public void removeModelListener(ModelListener<T> lst) { 
		listeners.remove(lst);
	}
	
	public Action createOkAction() { 
		return new AbstractAction("Ok") { 
			public void actionPerformed(ActionEvent e) { 
				ok();
			}
		};
	}
	
	public Action createCancelAction() { 
		return new AbstractAction("Cancel") { 
			public void actionPerformed(ActionEvent e) { 
				cancel();
			}
		};
	}
	
	public void ok() { 
		panel.saveToModel();
		T value = panel.getModelValue();
		dispose();
		for(ModelListener<T> lst : listeners) { 
			lst.modelChanged(value);
		}
		wakeWaiters();
	}
	
	private synchronized void wakeWaiters() { 
		notifyAll();
	}
	
	public synchronized T displayAndWait() { 
		display();
		try {
			wait();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		return panel.getModelValue();
	}
	
	public void cancel() { 
		dispose();
		wakeWaiters();
	}
	
	public void display() { 
		SwingUtilities.invokeLater(new Runnable() { 
			public void run() { 
				setLocation(100, 100);
				setVisible(true);
				pack();
			}
		});
	}   
}
