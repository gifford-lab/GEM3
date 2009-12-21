package edu.mit.csail.cgs.warpdrive.components;

import java.io.File;
import java.util.*;
import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.ewok.verbs.FastaWriter;

public class SaveRegionsAsFasta implements Runnable {
    
    private Collection<Region> regions;
    private JTextField both, up, down;

    public SaveRegionsAsFasta (Collection<Region> regions) {
        this.regions = regions;
        init();
    }
    public SaveRegionsAsFasta (Region r) {
        regions = new ArrayList<Region>();
        regions.add(r);
        init();
    }
    public SaveRegionsAsFasta(Iterator<Region> r) {
        regions = new ArrayList<Region>();
        while (r.hasNext()) {
            regions.add(r.next());
        }
        init();
    }
    public void init() {
        final JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(6,2));
        panel.add(new JLabel("Unstranded (eg Binding Events)"));
        panel.add(new JPanel());
        panel.add(new JLabel("bp to expand"));
        panel.add(both = new JTextField());
        panel.add(new JLabel("Stranded (eg Genes)"));
        panel.add(new JPanel());
        panel.add(new JLabel("bp to expand upstream"));
        panel.add(up = new JTextField());
        panel.add(new JLabel("bp to expand downstream"));
        panel.add(down = new JTextField());
        JButton ok, cancel;
        panel.add(ok = new JButton("OK"));
        panel.add(cancel = new JButton("Cancel"));
        final JFrame frame = new JFrame();
        final SaveRegionsAsFasta saver = this;
        both.setText("0");
        up.setText("0");
        down.setText("0");
        frame.getContentPane().add(panel);
        ok.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    Thread t = new Thread(saver);
                    t.start();
                    frame.dispose();
                }
            });
        cancel.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e){ 
                    frame.dispose();
                }
            });
        frame.pack();
        frame.setVisible(true);
    }
    public void run () {
        JFileChooser chooser;
        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
        int v = chooser.showSaveDialog(null);
        if(v == JFileChooser.APPROVE_OPTION) { 
            File f = chooser.getSelectedFile();
            int b = 0, u = 0, d = 0;
            try {
                b = Integer.parseInt(both.getText());
            } catch (NumberFormatException ex) {}
            try {
                u = Integer.parseInt(up.getText());
            } catch (NumberFormatException ex) {}
            try {
                d = Integer.parseInt(down.getText());
            } catch (NumberFormatException ex) {}
            try {
                FastaWriter fasta = new FastaWriter(f.getAbsolutePath());
                for (Region r : regions) {
                    if (r instanceof StrandedRegion) {
                        StrandedRegion sr = (StrandedRegion) r;
                        fasta.consume(sr.expand(u,d));
                    } else {
                        fasta.consume(r.expand(b,b));
                    }
                }
                fasta.finish();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
                
        }
    }
    

}
