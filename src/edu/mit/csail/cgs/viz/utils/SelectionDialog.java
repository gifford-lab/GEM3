/*
 * Created on Sep 29, 2005
 */
package edu.mit.csail.cgs.viz.utils;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;
import java.awt.event.*;

import java.util.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class SelectionDialog extends JDialog implements EventSource<SelectionEvent> {
    
    private Vector<String> options;
    private JButton ok, cancel;
    SelectionPanel<String> mainPanel;
    
    private EventSource.Default<SelectionEvent> src;
    
    void layoutDialog() {
        src = new EventSource.Default<SelectionEvent>();
        Container c = (Container)getContentPane();
        c.setLayout(new BorderLayout());
        
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1, 2, 5, 5));
        ok = new JButton("Ok");
        cancel = new JButton("Cancel");
        buttonPanel.add(ok);
        buttonPanel.add(cancel);
        
        ok.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) {
                SelectionEvent se = 
                    new SelectionEvent(SelectionDialog.this, 
                            SelectionEvent.OK, new Integer(getSelectionIndex()));
                SelectionDialog.this.dispose();
                src.fireEvent(se);                
            }
        });

        cancel.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) {
                SelectionEvent se = 
                    new SelectionEvent(SelectionDialog.this, SelectionEvent.CANCEL, null);
                SelectionDialog.this.dispose();
                src.fireEvent(se);
            }
        });

        mainPanel = new SelectionPanel<String>(options, 0);
        mainPanel.setBorder(new TitledBorder("Options:"));
        
        c.add(mainPanel, BorderLayout.CENTER);
        c.add(buttonPanel, BorderLayout.SOUTH);
        
        addWindowListener(new WindowAdapter() {
           public void windowClosed(WindowEvent e) {  
               SelectionEvent se = 
                   new SelectionEvent(SelectionDialog.this, SelectionEvent.CANCEL, null);
               src.fireEvent(se);               
           }
        });

        setVisible(true);
        pack();
        setLocation(getX() + 50, getY() + 50);
    }
    
    public int getSelectionIndex() { 
        return mainPanel.getSelectedIndex();
    }

    /**
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Collection<String> c) throws HeadlessException {
        super();
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Frame arg0, Collection<String> c) throws HeadlessException {
        super(arg0);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Frame arg0, boolean arg1, Collection<String> c) 
        throws HeadlessException {
        super(arg0, arg1);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Frame arg0, String arg1, Collection<String> c) 
    throws HeadlessException {
        super(arg0, arg1);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @param arg2
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Frame arg0, String arg1, boolean arg2, Collection<String> c)
            throws HeadlessException {
        super(arg0, arg1, arg2);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @param arg2
     * @param arg3
     */
    public SelectionDialog(Frame arg0, String arg1, boolean arg2,
            GraphicsConfiguration arg3, Collection<String> c) {
        super(arg0, arg1, arg2, arg3);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Dialog arg0, Collection<String> c) throws HeadlessException {
        super(arg0);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Dialog arg0, boolean arg1, Collection<String> c) 
        throws HeadlessException {
        super(arg0, arg1);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Dialog arg0, String arg1, Collection<String> c) 
        throws HeadlessException {
        super(arg0, arg1);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @param arg2
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Dialog arg0, String arg1, boolean arg2, Collection<String> c)
            throws HeadlessException {
        super(arg0, arg1, arg2);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /**
     * @param arg0
     * @param arg1
     * @param arg2
     * @param arg3
     * @throws java.awt.HeadlessException
     */
    public SelectionDialog(Dialog arg0, String arg1, boolean arg2,
            GraphicsConfiguration arg3, Collection<String> c) throws HeadlessException {
        super(arg0, arg1, arg2, arg3);
        options = new Vector<String>(c);
        layoutDialog();
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.EventSource#addEventListener(edu.mit.csail.cgs.utils.Listener)
     */
    public void addEventListener(Listener<SelectionEvent> el) {
        src.addEventListener(el);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.EventSource#removeEventListener(edu.mit.csail.cgs.utils.Listener)
     */
    public void removeEventListener(Listener<SelectionEvent> el) {
        src.removeEventListener(el);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.EventSource#hasListeners()
     */
    public boolean hasListeners() {
        return src.hasListeners();
    }

}
