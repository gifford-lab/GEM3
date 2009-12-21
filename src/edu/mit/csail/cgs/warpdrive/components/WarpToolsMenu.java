package edu.mit.csail.cgs.warpdrive.components;

import javax.swing.*;
import java.awt.event.*;
import edu.mit.csail.cgs.datasets.function.*;
import edu.mit.csail.cgs.tools.sequence.*;

public class WarpToolsMenu extends JMenu {

    private final RegionPanel panel;

    public WarpToolsMenu(RegionPanel p) {
        super("WarpDrive");
        panel = p;
        final RegionPanel thispanel = panel;
        JMenuItem item = new JMenuItem("Genome Browser");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    WarpOptionsFrame.main(new String[0]);
                }
            });
        item = new JMenuItem("Browse Weight Matrices");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    MotifDisplayPane.main(new String[0]);
                }
            });
        item = new JMenuItem("Gene Annotations");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    try {
                        DatabaseFunctionLoader loader = new DatabaseFunctionLoader();
                        GOAnnotationPanel panel = new GOAnnotationPanel(loader, thispanel.getRegion().getGenome());
                        JFrame f = new JFrame();
                        f.add(panel);
                        f.pack();
                        f.setVisible(true);
                    } catch (Exception ex) {
                        ex.printStackTrace();
                    }
                }
            });
        item = new JMenuItem("Sequence Browser");
        add(item);
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    if (panel == null) {
                        SequenceViewer seqview = new SequenceViewer();
                    }else {
                        SequenceViewer seqview = new SequenceViewer(panel.getRegion());
                    }
                }
            });
        if (panel != null) {
            item = new JMenuItem("Run CGH Analysis");
            add(item);
            item.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        if (panel != null) {
                            new CGHAnalysisPanel(panel);
                        } 
                    }
                }
                );
        }
    }


}
