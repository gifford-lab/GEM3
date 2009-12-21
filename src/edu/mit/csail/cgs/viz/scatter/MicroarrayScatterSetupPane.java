package edu.mit.csail.cgs.viz.scatter;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.util.Collection;
import java.util.ArrayList;
import java.sql.SQLException;
import edu.mit.csail.cgs.viz.components.ExptSelectPanel;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.datasets.species.*;

public class MicroarrayScatterSetupPane extends JTabbedPane implements ItemListener {

    private JPanel panelone, paneltwo;
    private JRadioButton cy5one, cy3one, ratioone, cy5two, cy3two, ratiotwo;
    private ButtonGroup radioone, radiotwo;
    private ExptSelectPanel exptselone, exptseltwo;
    private JCheckBox mvsa, logscale;
    private JComboBox species, genome;    
    private boolean handlingChange;
    
    public MicroarrayScatterSetupPane(Genome g) {
        handlingChange = false;
        exptselone = new ExptSelectPanel(g,
                                         false,
                                         true,
                                         false,
                                         false);
        exptseltwo = new ExptSelectPanel(g,
                                         false,
                                         true,
                                         false,
                                         false);
        cy5one = new JRadioButton("Cy5");
        cy3one = new JRadioButton("Cy3");
        ratioone = new JRadioButton("Ratio");
        cy5two = new JRadioButton("Cy5");
        cy3two = new JRadioButton("Cy3");
        ratiotwo = new JRadioButton("Ratio");
        radioone = new ButtonGroup();
        radiotwo = new ButtonGroup();
        radioone.add(cy5one);
        radioone.add(cy3one);
        radioone.add(ratioone);
        cy5two.setSelected(true);
        cy3one.setSelected(true);
        radiotwo.add(cy5two);
        radiotwo.add(cy3two);
        radiotwo.add(ratiotwo);
        JPanel rpone = new JPanel(), rptwo = new JPanel();
        rpone.add(cy5one);
        rpone.add(cy3one);
        rpone.add(ratioone);
        rptwo.add(cy5two);
        rptwo.add(cy3two);
        rptwo.add(ratiotwo);
        panelone = new JPanel();
        panelone.setLayout(new BorderLayout());
        paneltwo = new JPanel();
        paneltwo.setLayout(new BorderLayout());
        panelone.add(exptselone, BorderLayout.CENTER);
        paneltwo.add(exptseltwo, BorderLayout.CENTER);
        panelone.add(rpone, BorderLayout.SOUTH);
        paneltwo.add(rptwo, BorderLayout.SOUTH);
        
        JPanel optsPanel = new JPanel();
        species = new JComboBox();
        genome = new JComboBox();
        Collection<String> organisms = Organism.getOrganismNames();
        for (String o : organisms) {
            species.addItem(o);
        }
        if (g != null) {
            species.setSelectedItem(g.getSpecies());
            updateGenomeSelection();
            genome.setSelectedItem(g.getVersion());
        } else {
            species.setSelectedIndex(0);
        }
        species.addItemListener(this);
        genome.addItemListener(this);        
        optsPanel.setLayout(new GridLayout(3,2));
        optsPanel.add(new JLabel("Species"));
        optsPanel.add(species);
        optsPanel.add(new JLabel("Genome"));
        optsPanel.add(genome);
        JPanel junk = new JPanel();
        junk.setLayout(new BorderLayout());
        junk.add(optsPanel, BorderLayout.NORTH);
        junk.add(new JPanel(), BorderLayout.CENTER);
        
        mvsa = new JCheckBox("M vs A");
        logscale = new JCheckBox("Use Log Scale");
        optsPanel.add(mvsa);
        optsPanel.add(logscale);

        addTab("Options", junk);
        addTab("Channel One (X Axis)", panelone);
        addTab("Channel Two (Y Axis)", paneltwo);
    }
    public void setMvsA(boolean b) {mvsa.setSelected(b);}
    public void setLogScale(boolean b) {logscale.setSelected(b);}
    public void setExptOne(ExptLocator l) {
        exptselone.addToSelected(l);
    }
    public void setExptTwo(ExptLocator l) {
        exptseltwo.addToSelected(l);
    }
    
    private void updateGenomeSelection () {
        try {
            Organism org = new Organism(species.getSelectedItem().toString());
            Collection<String> genomes = org.getGenomeNames();
            genome.removeAllItems();
            for (String o : genomes) {
                genome.addItem(o);
            }
            genome.setSelectedIndex(0);                
        } catch (NotFoundException ex) {
            System.err.println("Couldn't find species " + species.getSelectedItem());
            ex.printStackTrace();
        }
    }
    public void itemStateChanged(ItemEvent e) {
        if (handlingChange) {return;}
        Object source = e.getItemSelectable();
        if (source == species) {
            updateGenomeSelection();
        }
        if (source == genome ||
            source == species) {
            synchronized(this) {
                if (!handlingChange) {
                    handlingChange = true;
                    updateExptSelection();
                    handlingChange = false;
                }
            }
        }        
    }
    private void updateExptSelection() {
        try {
            Organism org = new Organism(species.getSelectedItem().toString());
            Genome g = org.getGenome(genome.getSelectedItem().toString());
            exptselone.setGenome(g);
            exptseltwo.setGenome(g);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    /* Uses the selected experiments to create
       the dataset.  If multiple experiments
       are selected, then just takes an arbitrary one.  You can't
       really do anything more intelligen because the ExptSelectPanel
       returns a Collection that may not be ordered.
    */
    public Dataset2D parse() {
        Collection<ExptLocator> c1 = exptselone.getObjects();
        Collection<ExptLocator> c2 = exptseltwo.getObjects();
        try {
            ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
            int exptidone = -1, exptidtwo = -1;
            String labelone = "", labeltwo = "";
            System.err.println("Getting exptidone");
            for (ExptLocator l : c1) {
                ExptNameVersion env = (ExptNameVersion) l.getNameVersion();
                Experiment expt = loader.loadExperiment(env.getName(),
                                                        env.getVersion(),
                                                        env.getReplicate());
                exptidone = expt.getDBID();
                labelone = env.toString();
            }
            System.err.println(labelone + "\t" + exptidone);
            for (ExptLocator l : c2) {
                ExptNameVersion env = (ExptNameVersion) l.getNameVersion();
                Experiment expt = loader.loadExperiment(env.getName(),
                                                        env.getVersion(),
                                                        env.getReplicate());
                exptidtwo = expt.getDBID();
                labeltwo = env.toString();
            }
            System.err.println(labeltwo + "\t" + exptidtwo);

            PairedSQLData data = new PairedSQLData(exptidone, exptidtwo);
            int typeone = PairedSQLData.CHANNELONE, typetwo = PairedSQLData.CHANNELTWO;
            if (mvsa.isSelected()) {
                cy5one.setSelected(true);
                cy3two.setSelected(true);
            }
            if (cy5one.isSelected()) {
                typeone = PairedSQLData.CHANNELONE;
            } else if (cy3one.isSelected()) {
                typeone = PairedSQLData.CHANNELTWO;
            } else if (ratioone.isSelected()) {
                typeone = PairedSQLData.RATIO;
            }
            if (cy5two.isSelected()) {
                typetwo = PairedSQLData.CHANNELONE;
            } else if (cy3two.isSelected()) {
                typetwo = PairedSQLData.CHANNELTWO;
            } else if (ratiotwo.isSelected()) {
                typetwo = PairedSQLData.RATIO;
            }
            System.err.println("Typeone is " + typeone + " and tyypetwo is " + typetwo);
            float[][] array = data.getDataLimited( typeone, typetwo, 55000);
            if (mvsa.isSelected()) {
                for (int i = 0; i < array[0].length; i++) {
                    double a = array[0][i];
                    double b = array[1][i];
                    array[0][i] = (float)(.5 * Math.log(a*b));
                    array[1][i] = (float) Math.log(b/a);
                }
                String lone = labelone, ltwo = labeltwo;
                labelone = String.format(".5 * log(%s * %s", lone, ltwo);
                labeltwo = String.format("log(%s/%s)", ltwo, lone);
            } else if (logscale.isSelected()) {
                for (int i = 0; i < array[0].length; i++) {
                    array[0][i] = (float) Math.log(array[0][i]);
                    array[1][i] = (float) Math.log(array[1][i]);
                }
                labelone = "log(" +labelone + " " + PairedSQLData.columns[typeone] + ")";
                labeltwo = "log(" + labeltwo + " " + PairedSQLData.columns[typetwo] + ")";
            } else {
                labelone = labelone + " " + PairedSQLData.columns[typeone];
                labeltwo = labeltwo + " " + PairedSQLData.columns[typetwo];
            }
            System.err.println("labelone is " + labelone);
            System.err.println("labeltwo is " + labeltwo);
            return new Dataset2D(array,
                                 labelone,
                                 labeltwo);
                                 
        } catch (SQLException e) {
            throw new DatabaseException(e.toString(), e);
        } catch (NotFoundException e) {
            throw new DatabaseException(e.toString(), e);
        }              
    }       

}