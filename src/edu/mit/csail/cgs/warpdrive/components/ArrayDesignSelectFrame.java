package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.Collection;
import java.util.Iterator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.ewok.verbs.TiledRegionGenerator;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;

public class ArrayDesignSelectFrame extends JFrame implements ActionListener {

    private RegionList list;
    private JButton ok, cancel;    
    private JComboBox box;
    private DefaultComboBoxModel model;
    private Genome genome;

    public ArrayDesignSelectFrame(Genome g, RegionList list) {
        super();
        this.list = list;
        this.genome = g;
        try {
            model = new DefaultComboBoxModel();
            box = new JComboBox(model);        
            ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
            Collection<ArrayDesign> designs = loader.loadAllArrayDesigns(new Organism(g.getSpecies()));
            for (ArrayDesign d : designs) {
                model.addElement(d);
            }
            JPanel buttonPanel = new JPanel();
            ok = new JButton("OK");
            ok.addActionListener(this);
            cancel = new JButton("Cancel");
            cancel.addActionListener(this);
            buttonPanel.add(ok);
            buttonPanel.add(cancel);
            JPanel main = new JPanel();
            main.setLayout(new BorderLayout());
            main.add(box,BorderLayout.CENTER);
            main.add(buttonPanel,BorderLayout.SOUTH);
            add(main);
            pack();
            setVisible(true);
        } catch (Exception e) {
            e.printStackTrace();
        }        
    }
    public void actionPerformed (ActionEvent e) {
        if (e.getSource() == ok) {
            System.err.println("Creating populator");
            Populator p = new Populator(((ArrayDesign)model.getSelectedItem()).getName());
            System.err.println("Creating thread");
            Thread t = new Thread(p);
            System.err.println("Running thread");
            t.start();
            System.err.println("Disposing of this");
            this.dispose();
        } else if (e.getSource() == cancel) {
            this.dispose();
        }
    }

    private class Populator implements Runnable {
        private String designname;
        public Populator(String n) {
            designname = n;
        }
        public void run () {
            try {
                ChromosomeGenerator chromgen = new ChromosomeGenerator();
                TiledRegionGenerator tiled = new TiledRegionGenerator(designname,
                                                                      3000,
                                                                      3);
                JFrame frame = new JFrame();
                JPanel panel = new JPanel();
                JLabel label = new JLabel("Adding Tiled Regions...");
                panel.add(label);
                frame.add(panel);
                frame.pack();
                frame.setVisible(true);
                Iterator<Region> iter = new ExpanderIterator(tiled,chromgen.execute(genome));
                int added = 0;
                while (iter.hasNext()) {
                    list.addRegion(iter.next());
                    if (added++ % 50 == 0) {
                        label.setText("Adding Tiled Regions..." + added);
                    }
                }    
                frame.dispose();
            } catch (NotFoundException e) {
                e.printStackTrace();
            }            

        }


    }


}
