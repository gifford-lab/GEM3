package edu.mit.csail.cgs.viz.paintable;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.LinkedList;

import javax.swing.Action;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

public class PaintableFrame 
extends JFrame { 
    
    /**
     * Comment for <code>serialVersionUID</code>
     */
    private static final long serialVersionUID = 1L;
    
    private PaintablePanel fPanel;
    
    public PaintableFrame(String name, Paintable p) { 
        super(name);
        fPanel = new PaintablePanel(p);
        
        Container c = (Container)getContentPane();
        c.setLayout(new BorderLayout());
        c.add(fPanel, BorderLayout.CENTER);
        fPanel.setPreferredSize(new Dimension(1000, 500));
        
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setJMenuBar(createDefaultJMenuBar());
        setVisible(true);
        pack();
        setLocation(100, 100);
    }
    
    private JMenuBar createDefaultJMenuBar() { 
        JMenuBar jmb = new JMenuBar();
        JMenu menu; JMenuItem item;
        
        jmb.add((menu = new JMenu("File")));
        menu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) { 
                System.exit(0);
            }
        });
        
        
        LinkedList<Action> actions = collectActions();
        jmb.add((menu = new JMenu("Actions")));
        //menu.add(new JMenuItem(fPanel.getPaintable().getSaveImageAction()));
        for(Action a : actions) { 
            menu.add(new JMenuItem(a));
        }
        
        return jmb;
    }

	public PaintablePanel getPanel() { return fPanel; }
    
    private LinkedList<Action> collectActions() { 
        LinkedList<Action> lst = new LinkedList<Action>();
        Collection<Paintable> pCollect = fPanel.getPaintables();
        for(Paintable p : pCollect) { 
            lst.addAll(p.getPaintableActions());
        }

        return lst;
    }
}
