/*
 * Created on Sep 5, 2005
 */
package edu.mit.csail.cgs.utils.preferences;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class SelectionPreferencesPanel<X> extends PreferencesPanel {
    
    private JComboBox box;
    private DefaultComboBoxModel model;
    private Vector<IndexedPreferences<X>> indexed;
    private JButton prefsButton;

    public SelectionPreferencesPanel(int selected, Vector<Preferences<X>> prefs) {
        super();
        model = new DefaultComboBoxModel();
        indexed = new Vector<IndexedPreferences<X>>();
        for(int i = 0; i < prefs.size(); i++) { 
            IndexedPreferences<X> ip = new IndexedPreferences<X>(i, prefs.get(i));
            indexed.add(ip);
            model.addElement(ip);
        }
        model.setSelectedItem(indexed.get(selected));
        box = new JComboBox(model);
        prefsButton = new JButton("Preferences");
        
        prefsButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                IndexedPreferences<X> ip = 
                    (IndexedPreferences<X>)model.getSelectedItem();
                PreferencesPanel pp = ip.prefs.createPanel();
                PreferencesFrame pf = new PreferencesFrame(pp);
                pf.addEventListener(new FrameListener<X>(ip.prefs));
            }
        });
        
        setLayout(new BorderLayout());
        add(box, BorderLayout.NORTH);
        add(prefsButton, BorderLayout.SOUTH);
    }

    public void saveValues() { 
        super.saveValues();
        IndexedPreferences<X> sel = (IndexedPreferences<X>)model.getSelectedItem();
        for(int i = 0; i < indexed.size(); i++) { 
            if(sel.equals(indexed.get(i))) { 
                values.put("selected", indexed.get(i));
                return;
            }
        }
    }
}

class FrameListener<Y> implements Listener<PreferencesEvent> {
    
    private Preferences<Y> basePrefs;
    
    public FrameListener(Preferences<Y> p) { 
        basePrefs = p;
    }
    
    public void eventRegistered(PreferencesEvent pe) { 
        if(pe.getType() == PreferencesEvent.OK) { 
            PreferencesFrame pf = (PreferencesFrame)pe.getSource();
            PreferencesPanel pp = pf.getPanel();
            basePrefs.saveFromPanel(pp);
        }
    }
}

class IndexedPreferences<Y> { 
    public int index;
    public Preferences<Y> prefs;
    
    public IndexedPreferences(int i, Preferences<Y> p) { 
        index = i;
        prefs = p;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof IndexedPreferences)) { return false; }
        return index == ((IndexedPreferences)o).index;
    }
    
    public int hashCode() { 
        int code = 17;
        code += index; code *= 37;
        return code;
    }
    
    public String toString() { return prefs.toString(); }
}
