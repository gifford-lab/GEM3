package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.cgh.*;
import edu.mit.csail.cgs.viz.components.*;

public class CGHAnalysisPanel extends JPanel implements Runnable, ActionListener {

    private ExptSelectPanel exptpanel;
    private JPanel paramspanel, buttonpanel, bottompanel;
    private RegionPanel regionpanel;
    private Genome genome;

    private JTextField lowmean, lowstd, midmean, midstd, highmean, highstd, probemaxcount;    
    JButton ok, cancel;
    private JFrame paramsframe;

    public CGHAnalysisPanel (RegionPanel rp) {
        super();
        regionpanel = rp;
        genome = rp.getGenome();
        exptpanel = new ExptSelectPanel(genome);

        exptpanel.setPreferredSize(new Dimension(400,600));        

        lowmean = new JTextField(String.format("%f", CGHCallExpander.defaultParams[0][0]));
        lowstd = new JTextField(String.format("%f", Math.sqrt(CGHCallExpander.defaultParams[0][1])));

        midmean = new JTextField(String.format("%f", CGHCallExpander.defaultParams[1][0]));
        midstd = new JTextField(String.format("%f", Math.sqrt(CGHCallExpander.defaultParams[1][1])));

        highmean = new JTextField(String.format("%f", CGHCallExpander.defaultParams[2][0]));
        highstd = new JTextField(String.format("%f", Math.sqrt(CGHCallExpander.defaultParams[2][1])));
        
        paramspanel = new JPanel();
        paramspanel.setLayout(new GridLayout(5,3));
        paramspanel.add(new JLabel("State"));
        paramspanel.add(new JLabel("mean of log-ratio"));
        paramspanel.add(new JLabel("stddev of log-ratio"));
        paramspanel.add(new JLabel("Low"));
        paramspanel.add(lowmean);
        paramspanel.add(lowstd);

        paramspanel.add(new JLabel("Mid"));
        paramspanel.add(midmean);
        paramspanel.add(midstd);

        paramspanel.add(new JLabel("High"));
        paramspanel.add(highmean);
        paramspanel.add(highstd);

        paramspanel.add(new JLabel("Probe max count"));
        probemaxcount = new JTextField("-1");
        paramspanel.add(probemaxcount);

        buttonpanel = new JPanel();
        buttonpanel.setLayout(new GridLayout(1,2));
        ok = new JButton("OK");
        ok.addActionListener(this);
        buttonpanel.add(ok);
        cancel = new JButton("Cancel");
        cancel.addActionListener(this);
        buttonpanel.add(cancel);

        bottompanel = new JPanel();
        bottompanel.setLayout(new BorderLayout());
        bottompanel.add(paramspanel,BorderLayout.CENTER);
        bottompanel.add(buttonpanel,BorderLayout.SOUTH);

        setLayout(new BorderLayout());
        add(exptpanel, BorderLayout.CENTER);
        add(bottompanel, BorderLayout.SOUTH);

        paramsframe = new JFrame("CGH Parameters and Experiment Selection");
        paramsframe.getContentPane().add(this);
        paramsframe.setVisible(true);
        paramsframe.pack();
    }

    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == ok) {
            Thread t = new Thread(this);
            t.start();
            paramsframe.dispose();
        } else if (e.getSource() == cancel) {
            paramsframe.dispose();
        }
    }

    public void run() {
        Collection<ExptLocator> expts = exptpanel.getSelected();
        exptpanel.close();
        try {
            ChipChipDataset dataset = new ChipChipDataset(genome);
            ChromosomeGenerator chromgen = new ChromosomeGenerator();
            for (ExptLocator expt : expts) {
                ArrayList<ScoredRegion> output = new ArrayList<ScoredRegion>();
                ChipChipData data;
                try {
                    data = dataset.getData((ExptNameVersion)expt);
                } catch (ClassCastException e) {
                    continue;
                }

                try {
                    ((SQLData)data).setMaxCount(Integer.parseInt(probemaxcount.getText()));
                } catch (ClassCastException e) {
                    System.err.println("Couldn't set max probe count in CGH analysis");
                }

                CGHCallExpander cgh = new CGHCallExpander(data);
                cgh.setStateParameters(0,
                                       Double.parseDouble(lowmean.getText()),
                                       Double.parseDouble(lowstd.getText()));

                cgh.setStateParameters(1,
                                       Double.parseDouble(midmean.getText()),
                                       Double.parseDouble(midstd.getText()));

                cgh.setStateParameters(2,
                                       Double.parseDouble(highmean.getText()),
                                       Double.parseDouble(highstd.getText()));

                Iterator<Region> chroms = chromgen.execute(genome);
                while (chroms.hasNext()) {
                    Iterator<ScoredRegion> calls = cgh.execute(chroms.next());
                    while (calls.hasNext()) {
                        output.add(calls.next());
                    }
                }
                RegionListPanel rlp = new RegionListPanel(regionpanel, output);
                RegionListPanel.makeFrame(rlp,"CGH Results for " + expt.toString());
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }



    }





}