package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
import java.util.*;
import java.io.File;
import java.io.IOException;

import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.viz.DynamicAttribute;
import edu.mit.csail.cgs.warpdrive.*;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;

/* RegionFrame holds a RegionPanel
 */

public class RegionFrame extends JFrame {
    
    public static void main(String args[]) throws NotFoundException {
        WarpOptions opts = WarpOptions.parseCL(args);
        RegionFrame frame = new RegionFrame(opts);
    }

    public static Map<String,RegionFrame> genomeFrameMap;
    
    static { 
        genomeFrameMap = new HashMap<String,RegionFrame>();
    }
    
    public static void registerRegionFrame(RegionFrame f) { 
        String genome = f.panel.getGenome().getVersion();
        if(!genomeFrameMap.containsKey(genome)) { 
            genomeFrameMap.put(genome, f);
        } else { 
            System.err.println("There seems to be a duplicate RegionFrame \"" + genome + "\"");
        }
    }
    
    public static void unregisterRegionFrame(RegionFrame f) { 
        String genome = f.panel.getGenome().getVersion();
        if(!genomeFrameMap.containsKey(genome)) { throw new IllegalArgumentException(genome); }
        genomeFrameMap.remove(genome);        
    }
    
    public static RegionFrame findOrCreateRegionFrame(String genomeName) { 
        if(genomeFrameMap.containsKey(genomeName)) { 
            return genomeFrameMap.get(genomeName);
        } else { 
            WarpOptions opts = new WarpOptions(genomeName);
            RegionFrame f = new RegionFrame(opts);
            return f;
        }
    }
    
    private RegionPanel panel;
    private boolean imageraster;
    private int imageheight = 1200, imagewidth = 1600;

    public RegionFrame(WarpOptions opts) {
        setTitle(opts.species + " " + opts.genome);
        panel = new RegionPanel(opts);
        //getContentPane().add(new ImageCachingPanel(panel));
        getContentPane().add(panel);
        setJMenuBar(createDefaultJMenuBar(opts));
        setSize(600,400);
        setLocation(50,50);
        registerRegionFrame(this);
        setVisible(true);
        imageraster = true;

        this.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent arg0) {
                System.out.println("WindowClosing: " + arg0.toString());
                unregisterRegionFrame(RegionFrame.this);
                panel.handleWindowClosing();
            }
        });
    }
    
    public RegionPanel getRegionPanel() { return panel; }
    
    private JMenuBar createDefaultJMenuBar(WarpOptions opts) { 
        JMenuBar jmb = new JMenuBar();
        JMenu filemenu, imagemenu, navigationmenu, displaymenu, toolsmenu; 
        JMenuItem item;
        final RegionFrame thisframe = this;
        final RegionPanel thispanel = panel;
        jmb.add((filemenu = new JMenu("File")));
        filemenu.add((item = new JMenuItem("Configure Tracks")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {                    
                    WarpOptions current = panel.getCurrentOptions();
                    WarpOptionsFrame frame;
                    try {
                        frame = new WarpOptionsFrame(panel.getCurrentOptions());
                        frame.addPainterContainer(panel);
                    } catch (NotFoundException e1) {
                        e1.printStackTrace();
                    }
                }
            });
        filemenu.add((item = new JMenuItem("Add a track from a file")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {                    
                    thispanel.addTrackFromFile();
                }
            });
        filemenu.add((item = new JMenuItem("Save Region as FASTA")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new SaveRegionsAsFasta(thispanel.getRegion());
                }
            });
        filemenu.addSeparator();
        filemenu.add((item = new JMenuItem("Close")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    thisframe.dispose();
                }
            });
        filemenu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) { 
                    thispanel.close();
                    System.exit(0);
                }
            });        

        jmb.add((imagemenu = new JMenu("Image")));
        imagemenu.add((item = new JMenuItem("Save Image")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showSaveDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        boolean raster = imageraster;
                        if (f.getAbsolutePath().matches(".svg")) {
                            raster = false;
                        }
                        try {
                            panel.saveImage(f,imagewidth,imageheight,raster);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                    }
                } 
            });
        imagemenu.add((item = new JMenuItem("Image Settings")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new ImageConfigurationFrame(thisframe);
                }
            });
        jmb.add((navigationmenu = new JMenu("Navigation")));
        final JCheckBoxMenuItem linkeditem;
        navigationmenu.add((linkeditem = new JCheckBoxMenuItem("Link to Alignments")));
        linkeditem.setSelected(false);
        linkeditem.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    if (linkeditem.isSelected()) {
                        SpeciesAlignFrame.addRegionPanel(panel);
                        panel.setRegion(panel.getRegion());
                    } else {
                        SpeciesAlignFrame.removeRegionPanel(panel);
                    }
                }
            });


        navigationmenu.add((item = new JMenuItem("ChipChip Binding Scan Annotation")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    BindingEventAnnotationPanel beap = new BindingEventAnnotationPanel(thispanel, new ArrayList<BindingEvent>());
                    JFrame f = new BindingEventAnnotationPanel.Frame(beap);
                    f.pack();
                    new BindingScanSelectFrame(thispanel.getGenome(),beap);
                }
            });

        navigationmenu.add((item = new JMenuItem("Binding Event List")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thispanel,null);
                    RegionListPanel.makeFrame(rlp);
                    new BindingScanSelectFrame(thispanel.getGenome(),rlp);
                }
            });               

        navigationmenu.add((item = new JMenuItem("Array Tiled Regions")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thispanel,null);
                    RegionListPanel.makeFrame(rlp);
                    new ArrayDesignSelectFrame(thispanel.getGenome(),rlp);
                }
            });

        navigationmenu.add((item = new JMenuItem("Array Tiled Regions and Genes")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionAnnotationPanel beap = new RegionAnnotationPanel(thispanel,null);
                    JFrame f = new RegionAnnotationPanel.Frame(beap);
                    f.pack();
                    new ArrayDesignSelectFrame(thispanel.getGenome(),beap);
                }
            });

        navigationmenu.add((item = new JMenuItem("Open Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showOpenDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        java.util.List<Region> regions = RegionPanel.readRegionsFromFile(panel.getGenome(),f.getAbsolutePath());
                        RegionListPanel p = new RegionListPanel(panel,
                                                                regions);
                        RegionListPanel.makeFrame(p);
                    }
                    
                }
            });
        navigationmenu.add((item = new JMenuItem("New Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel p = new RegionListPanel(panel,
                                                            new ArrayList<Region>());
                    RegionListPanel.makeFrame(p);
                }
            });
        ButtonGroup group = new ButtonGroup();
        jmb.add(displaymenu = new JMenu("Display"));
        displaymenu.add((item = new JRadioButtonMenuItem("Screen")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.SCREEN);
                    thisframe.repaint();
                }
            });
        item.setSelected(true);
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Display")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.DISPLAY);
                    thisframe.repaint();
                }
            });
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Print")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.PRINT);
                    thisframe.repaint();
                }
            });
        group.add(item);
        jmb.add(new WarpToolsMenu(panel));
        
        
        if(BLASTInterfaceFactory.defaultFactory != null &&
           BLASTInterfaceFactory.defaultFactory.size() > 0) { 
            JMenu blastMenu = new JMenu("BLAST Against"); 

            try {
                final Genome base = Organism.findGenome(opts.genome);
                for(String genome : BLASTInterfaceFactory.defaultFactory.getGenomes()) { 
                    JMenuItem bitem = new JMenuItem(genome);
                    blastMenu.add(bitem);

                    final Genome target = Organism.findGenome(genome);
                    bitem.addActionListener(new ActionListener() {
                            public void actionPerformed(ActionEvent e) {
                                RegionFrame targetFrame;
                                if(RegionFrame.genomeFrameMap.containsKey(target.getVersion())) { 
                                    targetFrame = RegionFrame.genomeFrameMap.get(target.getVersion());
                                } else { 
                                    WarpOptions opts = new WarpOptions(target.getVersion());                                                                        
                                    opts.chrom = "1";
                                    opts.start = 10000;
                                    opts.stop = 20000;
                                    targetFrame = new RegionFrame(opts);
                                }                                
                                RegionPanel targetPanel = targetFrame.getRegionPanel();
                    
                                /* we could create a BindingEventAnnotationPanel here instead ... */
                                RegionListPanel regionPanel = new RegionListPanel(targetPanel, null);
                                JFrame listFrame = RegionListPanel.makeFrame(regionPanel, "BLAST Hits");                                
                                RegionList regionList = regionPanel;

                                BlastRegionRunnable runner = 
                                    new BlastRegionRunnable(base, target, BLASTInterfaceFactory.defaultFactory.getInterface(target));                    
                                runner.startInThread(panel.getRegion(),regionList);
                                
                                Region startRegion = null;
                                if(regionList.regionListSize() > 0) {
                                    startRegion= regionList.regionAt(0);
                                } else { 
                                    java.util.List<String> chroms = target.getChromList();
                                    String firstChrom = chroms.get(0);
                                    int start = 1;
                                    int end = Math.min(10000, target.getChromLength(firstChrom));
                                    startRegion = new Region(target, firstChrom, start, end);
                                }
                                targetPanel.setRegion(startRegion);
                            }
                        });
                }

                jmb.add(blastMenu);
            } catch (NotFoundException e1) {
                e1.printStackTrace();
            }
        }

        return jmb;
    }

    class ImageConfigurationFrame extends JFrame implements ActionListener {
        private RegionFrame parent;
        private JCheckBox rasterbox;
        private JTextField widthfield, heightfield;
        private JButton okbutton, cancelbutton;

        public ImageConfigurationFrame(RegionFrame p) {
            parent = p;
            JLabel boxlabel = new JLabel("Configure Save-as-image");
            rasterbox = new JCheckBox("Raster Image?",parent.imageraster);
            JLabel widthlabel = new JLabel("Width");
            JLabel heightlabel = new JLabel("Height");
            widthfield = new JTextField(Integer.toString(parent.imagewidth));
            heightfield = new JTextField(Integer.toString(parent.imageheight));
            okbutton = new JButton("OK");
            cancelbutton = new JButton("Cancel");
            okbutton.addActionListener(this);
            cancelbutton.addActionListener(this);
            
            JPanel toppanel = new JPanel();
            toppanel.setLayout(new BorderLayout());
            
            JPanel buttonpanel = new JPanel();
            buttonpanel.add(okbutton);
            buttonpanel.add(cancelbutton);
            toppanel.add(buttonpanel,BorderLayout.SOUTH);

            JPanel infopanel = new JPanel();
            infopanel.setLayout(new BorderLayout());
            infopanel.add(rasterbox,BorderLayout.NORTH);
            JPanel textpanel = new JPanel();
            textpanel.setLayout(new GridLayout(2,2));
            textpanel.add(widthlabel);
            textpanel.add(widthfield);
            textpanel.add(heightlabel);
            textpanel.add(heightfield);
            infopanel.add(textpanel,BorderLayout.CENTER);
            
            toppanel.add(infopanel,BorderLayout.CENTER);

            getContentPane().add(toppanel);
            setMinimumSize(new Dimension(150,150));
            setSize(getPreferredSize());
            pack();
            setVisible(true);
        }

        public void actionPerformed (ActionEvent e) {
            if (e.getSource() == okbutton) {
                parent.imageraster = rasterbox.isSelected();
                try {
                    parent.imageheight = Integer.parseInt(heightfield.getText());
                } catch (NumberFormatException ex) {
                }
                try {
                    parent.imagewidth = Integer.parseInt(widthfield.getText());
                } catch (NumberFormatException ex) {
                }


                this.dispose();
            } else if (e.getSource() == cancelbutton) {
                this.dispose();
            }
        }     
   
    }
}
