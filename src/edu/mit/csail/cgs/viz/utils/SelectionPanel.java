/*
 * Created on Oct 11, 2005
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.viz.utils;

import java.util.Vector;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;


public class SelectionPanel<X> extends JPanel {
    
    private Vector<X> options;
    private Vector<JRadioButton> buttons;
    private ButtonGroup group;
    
    public SelectionPanel(Vector<X> opts, int sel) {
        super();
        options = opts;
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        buttons = new Vector<JRadioButton>();
        group = new ButtonGroup();

        for(X opt : options) { 
            JRadioButton rb = new JRadioButton(opt.toString());
            buttons.add(rb);
            group.add(rb);
            add(rb);
        }
        
        buttons.get(sel).setSelected(true);
    }
    
    public int getSelectedIndex() { 
        for(int i = 0; i < buttons.size(); i++) { 
            if(buttons.get(i).isSelected()) { return i; }
        }
        return -1;
    }
    
    public X getSelectedOption() { return options.get(getSelectedIndex()); }
}