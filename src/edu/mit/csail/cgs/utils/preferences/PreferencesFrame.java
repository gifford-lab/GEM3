/*
 * Created on Aug 22, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

import edu.mit.csail.cgs.utils.EventSource;
import edu.mit.csail.cgs.utils.Listener;

/**
 * @author tdanford
 */
public class PreferencesFrame extends JFrame { 
    
    private PreferencesPanel pp;
    private EventSource.Default<PreferencesEvent> src;
	private JButton ok, cancel;
    
    public PreferencesFrame(PreferencesPanel temp) {
        super("Preferences");
        pp = temp;
        src = new EventSource.Default<PreferencesEvent>(this);
        
        Container c = (Container)getRootPane();
        c.setLayout(new BorderLayout());
        c.add(pp, BorderLayout.CENTER);

		JPanel buttons = new JPanel(); buttons.setLayout(new GridLayout(1, 2));
		ok = new JButton("Ok");
		cancel = new JButton("Cancel");
		buttons.add(ok); buttons.add(cancel);
		c.add(buttons, BorderLayout.SOUTH);

		ok.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				PreferencesEvent pe = new PreferencesEvent(PreferencesFrame.this, PreferencesEvent.OK);
				pp.saveValues();
				src.fireEvent(pe);	
				dispose();
			}
		});

		cancel.addActionListener(new ActionListener() { 
			public void actionPerformed(ActionEvent e) { 
				PreferencesEvent pe = new PreferencesEvent(PreferencesFrame.this, PreferencesEvent.CANCEL);
				src.fireEvent(pe);	
				dispose();
			}
		});

		addWindowListener(new WindowAdapter() { 
			public void windowClosed(WindowEvent e) { 
				PreferencesEvent pe = new PreferencesEvent(PreferencesFrame.this, PreferencesEvent.CANCEL);
				src.fireEvent(pe);	
			}
		});

        setSize(400, 400);
		setLocation(100, 100);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        setVisible(true);
		pack();
    }
    
    public PreferencesPanel getPanel() { return pp; }
    public void addEventListener(Listener<PreferencesEvent> l) { src.addEventListener(l); }
    public void removeEventListener(Listener<PreferencesEvent> l) { src.removeEventListener(l); }
}

