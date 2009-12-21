package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.sql.SQLException;
import java.io.File;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.renderer.category.*; 
import org.jfree.chart.renderer.xy.*; 
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.warpdrive.WarpOptions;
import edu.mit.csail.cgs.viz.paintable.*;
import edu.mit.csail.cgs.viz.charting.*;

/* Usage:
   java edu.mit.csail.cgs.warpdrive.components.ChipChipScatterFrame --species "$SC;Sigmav6" 
     --agilent "Sc WT poly-A RNA vs Cti6 poly-A RNA;no norm;3/24/08 (sigma 2-array set)" 
     [--mvsa]
     [--save foo.png]

*/

public class ChipChipScatterFrame extends JFrame {

    JFreeChart chart;
    Color color;

    public ChipChipScatterFrame(SQLData data,
                                String key,
                                Genome genome,
                                double sample,
                                Set<String> flags) throws SQLException {
        final ChipChipScatterFrame thisFrame = this;
        JMenuBar menubar = new JMenuBar();
        JMenu filemenu, colormenu;
        filemenu = new JMenu("File");
        JMenuItem item;
        filemenu.add((item = new JMenuItem("Close")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    thisFrame.dispose();
                }});
        filemenu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    System.exit(0);
                }});
        menubar.add(filemenu);

        colormenu = new JMenu("Color");
        colormenu.add((item = new JMenuItem("Set Color")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    thisFrame.chooseColor();
                }});
        menubar.add(colormenu);

        setJMenuBar(menubar);

        DefaultXYDataset dataset = new DefaultXYDataset();
        float[][] rawvalues = data.getRawValues();
        int tokeep = 0;
        boolean mvsa = flags.contains("mvsa");
        boolean log = flags.contains("log") && !mvsa;
        for (int i = 0; i < rawvalues[0].length; i++) {
            if (mvsa) {
                double ip =  rawvalues[0][i];
                double wce = rawvalues[1][i];
                rawvalues[1][i] = (float)Math.log(ip / wce);
                rawvalues[0][i] = (float)(.5 * Math.log(ip * wce));
            }
            if (log) {
                rawvalues[0][i] = (float)Math.log(rawvalues[0][i]);
                rawvalues[1][i] = (float)Math.log(rawvalues[1][i]);
            }
            if (Double.isNaN(rawvalues[0][i]) ||
                Double.isNaN(rawvalues[1][i])) {
                continue;
            }
            tokeep++;
        }
        tokeep = (int)(tokeep * sample);
        double[][] values = new double[2][tokeep];
        int j = 0;
        for (int i = 0; i < rawvalues[0].length && j < tokeep; i++) {
            if (Double.isNaN(rawvalues[0][i]) ||
                Double.isNaN(rawvalues[1][i])) {
                continue;
            }
            if (Math.random() > (.9 * sample )) {
                continue;
            }
            /* IMPORTANT: note that we switch the order here.  That's intentional to
               make the axes work out right.
            */
            values[1][j] = rawvalues[0][i];
            values[0][j] = rawvalues[1][i];
            j++;
        }

        dataset.addSeries(key, values);
        String xlabel, ylabel;
        if (mvsa) {
            ylabel = "M = log(cy5/cy3)";
            xlabel = "A = .5 * log(cy5 * cy3)";
        } else if (log) {
            xlabel = "ln(Cy3)";
            ylabel = "ln(Cy5)";
        } else {
            xlabel = "Cy3";
            ylabel = "Cy5";
        }
        chart = ChartFactory.createScatterPlot(key,ylabel,xlabel,dataset,PlotOrientation.VERTICAL,false,false,false);
        chart.setAntiAlias(true);
        ChartPanel panel = new ChartPanel(chart);
        Plot plot= chart.getPlot();

        color = new Color((float)0, 
                          (float)0,
                          (float)1.0,
                          (float).1);

        XYItemRenderer renderer = ((XYPlot)plot).getRenderer();
        renderer.setSeriesPaint(0, color);
//         renderer.setBaseStroke(new BasicStroke((float).5));
//         renderer.setSeriesStroke(0,new BasicStroke((float).5));
        renderer.setSeriesShape(0,new Rectangle(1,1));

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(panel,BorderLayout.CENTER);
        setSize(500,500);        
    }

    public void chooseColor() {
        final ChipChipScatterFrame thisFrame = this;
        final JFrame frame = new JFrame("Plot Color");
        frame.getContentPane().setLayout(new BorderLayout());
        
        final JSlider slider = new JSlider(0,
                                           255,
                                           color.getAlpha());
        frame.getContentPane().add(slider,
                                   BorderLayout.NORTH);

        
        final JColorChooser chooser = new JColorChooser(color);
        frame.getContentPane().add(chooser,
                                   BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridBagLayout());
        Dimension buttonSize = new Dimension(40,25);
        JButton ok = new JButton("OK");
        JButton cancel = new JButton("Cancel");
        ok.setMaximumSize(buttonSize);
        cancel.setMaximumSize(buttonSize);
        ok.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    Color t = chooser.getColor();
                    color = new Color(t.getRed(),
                                      t.getGreen(),
                                      t.getBlue(),
                                      slider.getValue());
                    Plot plot= chart.getPlot();
                    XYItemRenderer renderer = ((XYPlot)plot).getRenderer();
                    renderer.setSeriesPaint(0, color);
                    frame.dispose();
                }});
        cancel.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    frame.dispose();
                }});
        buttonPanel.add(ok);
        buttonPanel.add(cancel);
        
        frame.getContentPane().add(buttonPanel,
                                   BorderLayout.SOUTH);
        frame.setSize(500,500);
        frame.setVisible(true);
    }

    public static void main(String args[]) throws Exception {
        WarpOptions options = WarpOptions.parseCL(args);
        Genome g = new Genome(options.species,options.genome);
        ChipChipDataset agilent = new ChipChipDataset(g);        
        Set<String> flags = Args.parseFlags(args);        
        for (ExptNameVersion env : options.agilentdata) {
            ChipChipData data = agilent.getData(env);
            if (data instanceof SQLData) {
                ChipChipScatterFrame f = new ChipChipScatterFrame((SQLData) data,
                                                                  env.toString(),                                                                  
                                                                  g,
                                                                  .1,
                                                                  flags);
                if (flags.contains("save")) {
                    String fname = env.toString() + ".png";
                    JFreeChartPaintableAdapter adapter = new JFreeChartPaintableAdapter(f.chart);
                    adapter.saveImage(new File(fname),
                                      1600,
                                      1200,
                                      true);  
                    
                } else {
                    f.setVisible(true);
                }
                System.err.println("Created for " + env.toString());
            }
        }
    }
}