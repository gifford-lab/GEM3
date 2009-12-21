package edu.mit.csail.cgs.warpdrive.components;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import edu.mit.csail.cgs.datasets.motifs.*;

public class MotifDisplayPane extends JSplitPane {

    private MotifSelectPanel selectPanel;
    private JTable table;
    private MotifDrawingTableModel model;

    public static void main(String args[]) {
        final JFrame frame = new JFrame();
        final MotifDisplayPane mdp = new MotifDisplayPane();
        frame.setContentPane(mdp);
        frame.setSize(800,800);
        frame.setLocation(50,50);
        JMenuBar jmb = new JMenuBar();
        JMenu filemenu = new JMenu("File");
        jmb.add(filemenu);
        JMenuItem item;
        filemenu.add (item = new JMenuItem("Close"));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    frame.dispose();
                }
            });
        JMenu imagemenu = new JMenu("Image");
        jmb.add(imagemenu);
        imagemenu.add(item = new JMenuItem("Save All"));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    mdp.saveAll();
                }
            });
        jmb.add(new WarpToolsMenu(null));
        frame.setJMenuBar(jmb);
        frame.pack();
        frame.setVisible(true);
    }

    public MotifDisplayPane() {
        super(JSplitPane.HORIZONTAL_SPLIT);
        setDividerLocation(.5);
        selectPanel = new MotifSelectPanel();
        selectPanel.filter();
        JPanel buttonPanel = new JPanel();
        JButton addButton = new JButton("Show Motifs");
        final MotifDisplayPane mdp = this;
        addButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    mdp.showMotifs();
                }
            });
        Dimension buttonSize = new Dimension(30,20);
        addButton.setMaximumSize(buttonSize);
        buttonPanel.setLayout(new GridBagLayout());
        buttonPanel.add(addButton);

        JPanel leftside = new JPanel();
        leftside.setLayout(new BorderLayout());
        leftside.add(selectPanel,BorderLayout.CENTER);
        leftside.add(buttonPanel,BorderLayout.SOUTH);

        model = new MotifDrawingTableModel();
        table = new JTable(model);
        table.setRowHeight(100);
        table.setDefaultRenderer(WeightMatrix.class,new MotifDrawingRenderer());
        table.addMouseListener(new MouseAdapter() {
                public void mouseClicked(MouseEvent e) {
                    if (e.getButton() == MouseEvent.BUTTON3) {
                        int row = table.rowAtPoint(e.getPoint());
                        WeightMatrix wm = model.getObject(row);
                        JFileChooser chooser;
                        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                        int v = chooser.showSaveDialog(null);
                        if(v == JFileChooser.APPROVE_OPTION) { 
                            try {
                                File f = chooser.getSelectedFile();
                                BufferedImage im = 
                                    new BufferedImage(800, 200, BufferedImage.TYPE_INT_RGB);
                                Graphics g = im.getGraphics();
                                Graphics2D g2 = (Graphics2D)g;
                                g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
                                WeightMatrixPainter wmp = new WeightMatrixPainter();
                                g2.setColor(Color.WHITE);
                                g2.fillRect(0,0,800,200);
                                wmp.paint(wm,g2,0,0,800,200);
                                ImageIO.write(im, "png", f);
                            } catch (IOException ex) {
                                ex.printStackTrace();
                            }
                        }
                    }
                }
            });
        JScrollPane drawingPanel = new JScrollPane(table);
        add(new JScrollPane(leftside));
        add(drawingPanel);
    }

    public void showMotifs() {
        model.clear();
        for (WeightMatrix m : selectPanel.getObjects()) {
            model.addObject(m);
        }
    }
    public void saveAll() {
        JFileChooser chooser;
        chooser = new JFileChooser(new File(System.getProperty("user.dir")));
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        int v = chooser.showSaveDialog(null);
        if(v == JFileChooser.APPROVE_OPTION) { 
            try {
                File directory = chooser.getSelectedFile();
                for (WeightMatrix wm : selectPanel.getObjects()) {
                    String name = wm.toString().replaceAll("\\W","_");
                    File outfile = new File(directory,name + ".png");
                    BufferedImage im = 
                        new BufferedImage(800, 200, BufferedImage.TYPE_INT_RGB);
                    Graphics g = im.getGraphics();
                    Graphics2D g2 = (Graphics2D)g;
                    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
                    WeightMatrixPainter wmp = new WeightMatrixPainter();
                    g2.setColor(Color.WHITE);
                    g2.fillRect(0,0,800,200);
                    wmp.paint(wm,g2,0,0,800,200);
                    ImageIO.write(im, "png", outfile);
                }                
            } catch (IOException ex) {
                ex.printStackTrace();
            }
                
        }
    }
}

class MotifDrawingTableModel extends MotifTableModel {
    
    public int getColumnCount() {
        //        return 2;
        return 1;
    }
    public Class getColumnClass(int i) {
        //        if (i == 0) {
        //            return String.class;
        //        } 
        //        if (i == 1) {
        return WeightMatrix.class;
        //        }
        //        return null;
    }
    public String getColumnName(int i) {
        return "";
    }
    public Object getValueAt(int row, int c) {
        //        if (c == 0) {
        //            return getWeightMatrix(row).toString();
        //        } else if (c == 1){
        return getObject(row);
        //        } else {
        //            return null;
        //        }
    }
}

class MotifDrawingRenderer implements TableCellRenderer {

    public Component getTableCellRendererComponent(JTable table,
                                                   Object value,
                                                   boolean isSelected,
                                                   boolean hasFocus,
                                                   int row,
                                                   int column) {
        System.err.println("Creating MotifCellRenderer for " + row + "," + column);
        return new MotifCellRenderer(table,value,row,column);
    }
}
class MotifCellRenderer extends JPanel {
    private JTable table;
    private Object value;
    private int row, column;
    private WeightMatrixPainter painter;
    public MotifCellRenderer (JTable table, Object value, int row, int column) {
        this.table = table;
        this.value = value;
        this.row = row;
        this.column = column;
        painter = new WeightMatrixPainter();
    }
    public void paintComponent(Graphics g) {
        if (value instanceof WeightMatrix) {
            Rectangle d = table.getCellRect(row,column,false);
            System.err.println("painting in " + d);
            painter.paint((WeightMatrix)value,
                          g,
//                           (int)d.getX(),
//                           (int)d.getY(),
//                           (int)(d.getX() + d.getWidth()),
//                           (int)(d.getY() + d.getHeight()));
                          0,0,(int)d.getWidth(),(int)d.getHeight());
        }
                          
    }
    public void paintComponent(Graphics g, int x, int y, int width, int height) {
        if (value instanceof WeightMatrix) {
            painter.paint((WeightMatrix)value,
                          g,
                          x,y,x+width,y+height);
        }
    }
}
