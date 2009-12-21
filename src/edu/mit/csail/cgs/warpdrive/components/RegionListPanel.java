package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import java.util.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;

public class RegionListPanel extends JPanel implements ActionListener, KeyListener, MouseListener, RegionList, Runnable {

    private RegionPanel regionpanel;
    private DefaultListModel listmodel;
    private JList regionlist;
    private JButton goButton, nextButton, prevButton, saveRegionButton;
    private JScrollPane posScrollPane;
    private Vector unadded;
    private boolean pendingadds;

    public RegionListPanel(RegionPanel panel, 
                           java.util.List<? extends Region> regions) {
        init();
        regionpanel = panel;
        if (regions != null) {
            for (int i = 0; i < regions.size(); i++) {
                listmodel.addElement(regions.get(i));
            }
        }
        unadded = new Vector();
        pendingadds = false;
    }

    public static JFrame makeFrame(RegionListPanel panel) {
        return makeFrame(panel,"Regions to View");
    }

    public static JFrame makeFrame(final RegionListPanel panel, String title) {
        final JFrame frame = new JFrame(title);
        frame.addKeyListener(panel);
        Container content = frame.getContentPane();
        content.setLayout(new BorderLayout());
        content.add(panel,BorderLayout.CENTER);
        content.add(new JLabel(title),BorderLayout.NORTH);
        JMenuBar jmb = new JMenuBar();
        JMenu menu = new JMenu("File");
        jmb.add(menu);
        JMenuItem item = new JMenuItem("Save");
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    panel.saveRegionList();
                }
            });
        menu.add(item);
        item = new JMenuItem("Save Tab-Separated Regions");
        item.addActionListener(new ActionListener() { 
            public void actionPerformed(ActionEvent e) { 
                panel.saveTabSeparatedRegionList();
            }
        });
        menu.add(item);
        item = new JMenuItem("Save as FASTA");
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    ArrayList<Region> regions = new ArrayList<Region>();
                    for (int i = 0; i < panel.listmodel.size(); i++) {
                        regions.add((Region)panel.listmodel.get(i));
                    }
                    new SaveRegionsAsFasta(regions);
                }
            });
        item = new JMenuItem("Close");
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    frame.dispose();
                }
            });
        menu.add(item);
        frame.setJMenuBar(jmb);
        frame.pack();
        frame.setVisible(true);
        return frame;
    }

    public void addRegion(Region r) {       
        synchronized(unadded) {
            boolean oldpending = pendingadds;
            unadded.add(r);
            pendingadds = true;
            if (!oldpending) {
                SwingUtilities.invokeLater(this);
            }
        }        
    }
    public void run() {
        synchronized(unadded) {
            for (Object r : unadded) {
                listmodel.addElement(r);
            }
            unadded.clear();
            pendingadds = false;
        }
    }
    public int regionListSize() {return listmodel.size();}
    public Region regionAt(int i) {return (Region)listmodel.getElementAt(i);}
    private void init() {
        listmodel = new DefaultListModel();
        regionlist = new JList(listmodel);
        regionlist.setLayoutOrientation(JList.VERTICAL);
        posScrollPane = new JScrollPane(regionlist);
        posScrollPane.setMinimumSize(new Dimension(150,300));
        JPanel buttonPanel = new JPanel();
        saveRegionButton = new JButton("Save Region to List");
        goButton = new JButton("Go");
        nextButton = new JButton("Next");
        prevButton = new JButton("Prev");
        buttonPanel.add(prevButton);
        buttonPanel.add(goButton);
        buttonPanel.add(nextButton);
        buttonPanel.add(new JPanel());
        buttonPanel.add(saveRegionButton);
        setLayout(new BorderLayout());
        add(posScrollPane,BorderLayout.CENTER);
        add(buttonPanel,BorderLayout.SOUTH);
        prevButton.addActionListener(this);
        goButton.addActionListener(this);
        nextButton.addActionListener(this);
        saveRegionButton.addActionListener(this);
        addKeyListener(this);
        regionlist.addMouseListener(this);
        regionlist.addKeyListener(this);
        this.setMinimumSize(new Dimension(200,300));
    }

    public void go() {
        Region r = (Region)regionlist.getSelectedValue();
        if (r == null) {return;}
        if (Math.abs(r.getStart() - r.getEnd()) < 200) {
            Region old = regionpanel.getRegion();
            int size = Math.abs(old.getStart() - old.getEnd());
            int center = (r.getStart() + r.getEnd())/2;
            r = new Region(r.getGenome(),
                           r.getChrom(),
                           center - size/2,
                           center + size/2);
        }
        regionpanel.setRegion(r);
    }

        public void first() {
        regionlist.ensureIndexIsVisible(0);
        regionlist.setSelectionInterval(0,0);
    }
    private void next() {
        move(1);
    }
    private void prev() {
        move(-1);
    }
    private void move (int offset) {
        int curpos = regionlist.getSelectedIndex();
        int newpos = curpos + offset;
        if ((curpos < 0) && (newpos < 0)) {
            return;
        }
        if (newpos < 0) {newpos = 0;}
        if (newpos >= listmodel.getSize()) {
            newpos = listmodel.getSize() - 1;
        }
        regionlist.ensureIndexIsVisible(newpos);
        regionlist.setSelectionInterval(newpos,newpos);
    }
    public void saveRegion() {
        addRegion(regionpanel.getRegion());
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == goButton) {
            go();            
        }else if (e.getSource() == nextButton) {
            next();
            go();
        }else if (e.getSource() == prevButton) {
            prev();
            go();
        } else if (e.getSource() == posScrollPane) {
            go();
        } else if (e.getSource() == saveRegionButton) {
            saveRegion();
        }
    }
    public void keyPressed(KeyEvent e) {
        System.err.println("KEY Typed " + e.getKeyCode());
        switch (e.getKeyCode()) {
        case KeyEvent.VK_N:
            next();
            break;
        case KeyEvent.VK_P:
            prev();
            break;
        case KeyEvent.VK_G:
            go();
            break;
        case KeyEvent.VK_ENTER:
            go();
            break;
        }        }
    public void keyReleased(KeyEvent e) {}
    public void keyTyped(KeyEvent e) {
    }
    public void mouseClicked(MouseEvent e) {
        if (e.getClickCount() == 2) {
            int newindex = regionlist.locationToIndex(e.getPoint());
            int currindex = regionlist.getSelectedIndex();
            move(newindex - currindex);
            go();
        }
    }
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}

    public void saveRegionList() {
        JFileChooser chooser;
        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
        int v = chooser.showSaveDialog(null);
        if(v == JFileChooser.APPROVE_OPTION) { 
            try {
                File f = chooser.getSelectedFile();
                PrintWriter writer = new PrintWriter(f);
                for (int i = 0; i < listmodel.getSize(); i++) {
                    writer.println(listmodel.get(i).toString());
                }
                writer.close();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }

    public void saveTabSeparatedRegionList() {
        JFileChooser chooser;
        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
        int v = chooser.showSaveDialog(null);
        if(v == JFileChooser.APPROVE_OPTION) { 
            try {
                File f = chooser.getSelectedFile();
                PrintWriter writer = new PrintWriter(f);
                for (int i = 0; i < listmodel.getSize(); i++) {
                    Region r = (Region)listmodel.get(i);
                    writer.print(r.getChrom() + "\t" + r.getStart() + "\t" + r.getEnd());
                    if(r instanceof BindingEvent) { 
                        BindingEvent b = (BindingEvent)r;
                        writer.print("\t" + b.getSize() + "\t" + b.getConf());
                    }
                    writer.println();
                }
                writer.close();
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }
}
