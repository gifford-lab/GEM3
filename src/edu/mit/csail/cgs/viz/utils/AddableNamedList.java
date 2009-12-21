/*
 * Created on Aug 21, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import edu.mit.csail.cgs.utils.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * @author tdanford
 */
public class AddableNamedList<X> extends NamedListPanel<X> implements ListPanelEventSource {

    private ObjectChooser<NamedValue<X>> chooser;
    private JButton addButton, deleteButton;
    private ListPanelEventSource.Default src;
    
    public AddableNamedList(ObjectChooser<NamedValue<X>> ch) {
        super();
        chooser = ch;
        src = new ListPanelEventSource.Default(this);
        addButton = new JButton("Add");
        deleteButton = new JButton("Delete");
        
        deleteButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                src.firePanelEvent(ListPanelEvent.REMOVE, this);
                removeSelectedValues();
            }
        });
        
        addButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                src.firePanelEvent(ListPanelEvent.ADD, this);
                if(chooser != null) { 
                    NamedValue<X> val = chooser.choose();
                    addNamedValue(val);
                }
            }
        });
        
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1, 2));
        buttonPanel.add(addButton);
        buttonPanel.add(deleteButton);
        
        add(buttonPanel, BorderLayout.SOUTH);
    }
    
    public void addListPanelListener(ListPanelListener l) { src.addListPanelListener(l); }
    public void removeListPanelListener(ListPanelListener l) { src.removeListPanelListener(l); }
}
