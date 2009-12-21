package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;

import edu.mit.csail.cgs.viz.colors.ColorSet;
import edu.mit.csail.cgs.warpdrive.paintable.SpeciesAlignPainter;
import edu.mit.csail.cgs.warpdrive.paintable.SpeciesAlignThinPainter;
import edu.mit.csail.cgs.warpdrive.model.SpeciesAlignModel;
import edu.mit.csail.cgs.utils.*;

public class SpeciesAlignFrame extends JFrame implements Listener<EventObject> {

    static SpeciesAlignFrame theFrame;
    private SpeciesAlignModel model;
    private SpeciesAlignPainter painter;
    private JMenu alignmentMenu;
    private JLabel alignmentLabel;
    private ArrayList<JCheckBoxMenuItem> menuitems;

    public static void addRegionPanel(RegionPanel p) {
        if (theFrame == null) {
            theFrame = new SpeciesAlignFrame();
        }
        theFrame.addPanel(p);        
    }
    
    public static void removeRegionPanel(RegionPanel p) {
        if (theFrame == null) {return;}
        theFrame.removePanel(p);
    }

    public SpeciesAlignFrame () {
        super();
        model = new SpeciesAlignModel();
        model.addEventListener(this);
        Thread t = new Thread(model);
        t.start();
        painter = new SpeciesAlignPainter(model);
        JPanel content = new JPanel();
        content.setLayout(new BorderLayout());
        alignmentLabel = new JLabel(model.getAlignment());
        alignmentLabel.setHorizontalAlignment( javax.swing.SwingConstants.CENTER);
        content.add(alignmentLabel,BorderLayout.NORTH);
        content.add(new SpeciesAlignPanel(painter),BorderLayout.CENTER);
        getContentPane().add(content);
        setTitle("alignment");
        setSize(300,200);
        setLocation(200,200);
        setVisible(true);
        JMenuBar jmb = new JMenuBar();
        alignmentMenu = new JMenu("Alignment Version");
        jmb.add(alignmentMenu);
        setJMenuBar(jmb);
        menuitems = new ArrayList<JCheckBoxMenuItem>();
    }

    public void addPanel(RegionPanel p) {
        model.addRegionPanel(p);
        SpeciesAlignThinPainter satp = new SpeciesAlignThinPainter(model,
                                                                   p.getGenome());
        satp.setLabel("Alignments");
        satp.addEventListener(p);
        model.addEventListener(satp);
        p.addPainter(satp);
        System.err.println("Added Thin painter to " + p);
        populateMenu();
    }    
    public void removePanel(RegionPanel p) {
        model.removeRegionPanel(p);
        p.removeTrack("Alignments");
        populateMenu();
    }
    
    class AlignmentActionListener implements ActionListener {
        private String alignment;
        public AlignmentActionListener (String a) {
            alignment = a;
        }
        public void actionPerformed(ActionEvent e) {
            model.setAlignment(alignment);
            alignmentLabel.setText(alignment);            
            for (JCheckBoxMenuItem item : menuitems) {
                if (item.getText().equals(alignment)) {
                    item.setSelected(true);
                } else {
                    item.setSelected(false);
                }
            }
        }
    }
    public void populateMenu() {
        alignmentMenu.removeAll();
        Set<String> methods = model.getAlignments();
        String current = model.getAlignment();
        for (String s : methods) {
            System.err.println(" ==> adding to AlignmentMenu " + s);
            JCheckBoxMenuItem item = new JCheckBoxMenuItem(s);
            menuitems.add(item);
            if (s.equals(current)) {
                item.setSelected(true);
            } else {
                item.setSelected(false);
            }
            item.addActionListener(new AlignmentActionListener(s));
            alignmentMenu.add(item);
        }
    }

    public void eventRegistered(EventObject e) {
        if (e.getSource() == model) {
            painter.setRegion(model.getRegion());
            repaint();
        }
    }

    class SpeciesAlignPanel extends JPanel implements MouseListener {
        private SpeciesAlignPainter painter;
        
        public SpeciesAlignPanel (SpeciesAlignPainter p) {
            super();
            addMouseListener(this);
            painter = p;
            p.setLabel("SpeciesAlignement");
        }
        public void paintComponent(Graphics g) {
            g.setColor(Color.WHITE);
            g.fillRect(getX(),getY(),getWidth(),getHeight());
            painter.paintItem((Graphics2D)g,
                              getX(),
                              getY(),
                              getX() + getWidth(),
                              getY() + getHeight());
                              
        }
         public void mouseClicked(MouseEvent e) {
//             if (e.getButton() == MouseEvent.BUTTON3) {
//                 painter.configurePaintable(this);
//             }
         }
        public void mouseEntered(MouseEvent e) {}

        public void mouseExited(MouseEvent e) {}

        public void mousePressed(MouseEvent e) {}
        
        public void mouseReleased(MouseEvent e) {}
    }

}

