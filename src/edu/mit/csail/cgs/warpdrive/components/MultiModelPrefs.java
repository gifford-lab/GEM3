package edu.mit.csail.cgs.warpdrive.components;

/** totally ripped off from edu.mit.csail.cgs.viz.eye.ModelPrefs but I want
 * to include multiple panels with different types without making Tim scream
 * when he sees what I did to his code...
 *
 * Also, this version doesn't do listeners since that's more complicated.
 */

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;
import java.lang.reflect.*;
import edu.mit.csail.cgs.viz.eye.*;
import edu.mit.csail.cgs.utils.models.*;

public class MultiModelPrefs extends JFrame {
	
	private ArrayList<PrefsPanel<Model>> panels;
	private JButton ok, cancel;
    private JPanel regionpanel;
	
	public MultiModelPrefs(Collection<? extends Model> models, JPanel rp) { 
		super("Preferences");
        regionpanel = rp;
        Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
        panels = new ArrayList<PrefsPanel<Model>>();
		
        JPanel mainpanel = new JPanel();
        GridBagLayout gridbag = new GridBagLayout();
        mainpanel.setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();        
        constraints.weightx = 1.0;
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        int height = 0;
        for (Model m : models) {
            if (m.getFields().isEmpty()) {
                continue;
            }
            PrefsPanel panel = new PrefsPanel<Model>(m);
            JLabel label = new JLabel(m.getClass().toString().replaceAll("^.*\\.",""));            
            height += (int)panel.getPreferredSize().getHeight() + 50;            
            gridbag.setConstraints(label,constraints);
            mainpanel.add(label);
            gridbag.setConstraints(panel,constraints);
            mainpanel.add(panel);
            panels.add(panel);
        }
        if (!panels.isEmpty()) {
            JScrollPane pane = new JScrollPane(mainpanel);
            pane.setPreferredSize(new Dimension(700,height));
//            System.err.println("preferred size is " + mainpanel.getPreferredSize());
            c.add(pane, BorderLayout.CENTER);
            JPanel buttons = new JPanel();
            buttons.setLayout(new FlowLayout());
            buttons.add(ok = new JButton(createOkAction()));
            buttons.add(cancel = new JButton(createCancelAction()));
            gridbag.setConstraints(buttons,constraints);
            c.add(buttons, BorderLayout.SOUTH);
            
            setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        } else {
            this.dispose();
        }
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
        for (PrefsPanel p : panels) {
            p.saveToModel();            
        }
        dispose();
		wakeWaiters();
        regionpanel.repaint();
	}
	
	private synchronized void wakeWaiters() { 
		notifyAll();
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


