package edu.mit.csail.cgs.warpdrive.components;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import java.sql.SQLException;
import java.util.Collection;
import java.util.ResourceBundle;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.*;

import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipDataset;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.datasets.general.NamedTypedRegion;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.viz.components.BindingScanSelectPanel;
import edu.mit.csail.cgs.viz.components.ExptSelectPanel;
import edu.mit.csail.cgs.warpdrive.WarpOptions;
import edu.mit.csail.cgs.ewok.RegionExpanderFactoryLoader;
import edu.mit.csail.cgs.ewok.PeakCallerFactoryLoader;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;

public class WarpOptionsPane 
	extends JTabbedPane 
	implements ItemListener, ActionListener, Closeable {
    // regexes get special handling at the moment because there's no gui component for them,
    // so cache them if neccessary
    private HashMap<String,String> regexes;
    private boolean handlingChange, closed;
    private RegionExpanderFactoryLoader<Gene> gfLoader;
    private RegionExpanderFactoryLoader<NamedTypedRegion> annotLoader;
    private PeakCallerFactoryLoader pcLoader;

    private JPanel speciesLocationPanel,
        annotationsPanel,
        chipChipPanel,
        chipSeqPanel,
        chipSeqAnalysisPanel,
        pairedChipSeqPanel,
        optionsPanel,
        peakPanel, 
        exprPanel;
    private WarpOptions createdFrom;
    private MotifScanSelectPanel motifScanPanel;
    private MotifSelectPanel motifPanel;
    
    // species/location tab
    private JComboBox species, genome;
    private JTextField position, gene;
    private JLabel specieslabel, genomelabel, positionlabel, genelabel;

    // options tab
    private JCheckBox relative, hash, common, seqletters, oldchipseq;
    
    // peak tab
    private BindingScanSelectPanel bindingSelect;
    
    // chipseq tab
    private ChipSeqSelectPanel chipSeqSelect, pairedChipSeqSelect;

    // annotations tab
    private JList genes, ncrnas, otherfeats;
    private DefaultListModel genesmodel, ncrnasmodel, otherfeatsmodel;
    private JCheckBox polyA, gccontent, pyrpurcontent, cpg, regexmatcher;
    private JLabel geneslabel, ncrnaslabel, otherfeatslabel;
    
    //private ExptTreeSelectPanel exptSelect;
    private ExptSelectPanel exptSelect;
    
    // expression tab
    private ExprExperimentSelectPanel exprSelect;

    // motif tab
    private JList motifs;

    // file-based tracks
    private FileBasedTracksPanel filetracks;

    public WarpOptionsPane () throws NotFoundException {
        super();
        gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
        pcLoader = new PeakCallerFactoryLoader();
        annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");
        handlingChange = true;
        init();
        handlingChange = false;
        closed = false;
        setSpeciesGenomeDefaults();
    }

    public WarpOptionsPane(String species, String genome) throws NotFoundException {
        super();
        gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
        pcLoader = new PeakCallerFactoryLoader();
        annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");
        handlingChange = true;
        init();
        handlingChange = false;
        closed = false;
        if (genome != null && species != null) {
            /* temporarily remove the item listener for species
               so that we don't trigger two updates for experiment selection
               and such.  Everything will be updated when we set the genome
            */
            this.genome.removeItemListener(this);
            this.species.removeItemListener(this);
            this.species.setSelectedItem(species);
            updateGenomeSelection();
            this.genome.setSelectedItem(genome);
            updateExptSelection();
            this.species.addItemListener(this);
            this.genome.addItemListener(this);
        } else if (species != null) {
            this.species.setSelectedItem(species);
        }

        regexes = null;
    }

    public WarpOptionsPane(WarpOptions opts) throws NotFoundException {
        super();
        gfLoader = new RegionExpanderFactoryLoader<Gene>("gene");
        pcLoader = new PeakCallerFactoryLoader();
        annotLoader = new RegionExpanderFactoryLoader<NamedTypedRegion>("annots");
        closed = false;
        init(opts);
        if (opts.genome == null) {
            setSpeciesGenomeDefaults();
        }        
        regexes = opts.regexes;
    }
    
    public boolean isClosed() { return closed; }
    
    public void close() { 
        exptSelect.close();
        exprSelect.close();
        chipSeqSelect.close();
        pairedChipSeqSelect.close();
    	bindingSelect.close(); 
    	closed = true; 
    }
    
    private Genome loadGenome() { 
        String orgName = (String)species.getSelectedItem();
        String genName = (String)genome.getSelectedItem();
        if(orgName != null && genName != null) { 
            Organism org;
            try {
                org = Organism.getOrganism(orgName);
                Genome gen = org.getGenome(genName);
                return gen;
            } catch (NotFoundException e) {
                e.printStackTrace();
                throw new IllegalArgumentException(orgName + ":" + genName);
            }
        }
        return null;
    }

    private void init() throws NotFoundException {
        speciesLocationPanel = new JPanel();
        annotationsPanel = new JPanel();
        chipChipPanel = new JPanel();       
        chipSeqPanel = new JPanel();
        chipSeqAnalysisPanel = new JPanel();
        pairedChipSeqPanel = new JPanel();
        optionsPanel = new JPanel();
        peakPanel = new JPanel();
        exprPanel = new JPanel();
        motifScanPanel = new MotifScanSelectPanel();
        motifPanel = new MotifSelectPanel();

        // First tab lets the user select the species, genome version,
        // and the genomic coordinates they want to view
        specieslabel = new JLabel("Species");
        genomelabel = new JLabel("Genome Version");
        positionlabel = new JLabel("Genome position\nto view");
        genelabel = new JLabel("Gene to view");
        species = new JComboBox();
        genome = new JComboBox();
        position = new JTextField();
        position.setText("1:10000-20000");
        gene = new JTextField();

        /* need to fill the species and genome boxes here */
        Collection<String> organisms = Organism.getOrganismNames();
        for (String o : organisms) {
            species.addItem(o);
        }

        species.setSelectedIndex(0);
        updateGenomeSelection();
        Organism org = new Organism(species.getSelectedItem().toString());
        Collection<String> genomes = org.getGenomeNames();
        genome.removeAllItems();
        for (String o : genomes) {
            genome.addItem(o);
        }
        Genome g = null;
        if(genome.getModel().getSize() > 0) { 
            genome.setSelectedIndex(0);
            String gname = (String)genome.getSelectedItem();
            g = Organism.findGenome(gname);
        } else {
            // umm, no genomes is bad and will break other stuff
        }        

        speciesLocationPanel.setLayout(new GridLayout(4,2));
        speciesLocationPanel.add(specieslabel);
        speciesLocationPanel.add(species);
        speciesLocationPanel.add(genomelabel);
        speciesLocationPanel.add(genome);
        speciesLocationPanel.add(positionlabel);
        speciesLocationPanel.add(position);
        speciesLocationPanel.add(genelabel);
        speciesLocationPanel.add(gene);
        species.addItemListener(this);
        genome.addItemListener(this);        
        
        // chipChip tab
        
        //exptSelect = new ExptTreeSelectPanel(null);
        exptSelect = new ExptSelectPanel(null);

        chipChipPanel.setLayout(new BorderLayout());
        chipChipPanel.add(exptSelect, BorderLayout.CENTER);
        
        // expression tab
        exprSelect = new ExprExperimentSelectPanel();
        
        exprPanel.setLayout(new BorderLayout());
        exprPanel.add(exprSelect, BorderLayout.CENTER);
        
        // chipseq tab
        chipSeqSelect = new ChipSeqSelectPanel();        
        chipSeqPanel.setLayout(new BorderLayout());
        chipSeqPanel.add(chipSeqSelect, BorderLayout.CENTER);

        pairedChipSeqSelect = new ChipSeqSelectPanel();
        pairedChipSeqPanel.setLayout(new BorderLayout());
        pairedChipSeqPanel.add(pairedChipSeqSelect, BorderLayout.CENTER);
        
        // chipseq analysis
        chipSeqAnalysisSelect = new ChipSeqAnalysisSelectPanel();
        chipSeqAnalysisPanel.setLayout(new BorderLayout());
        chipSeqAnalysisPanel.add(chipSeqAnalysisSelect, BorderLayout.CENTER);

        // peak tab
        peakPanel.setLayout(new BorderLayout());
        try {
            bindingSelect = new BindingScanSelectPanel();
            peakPanel.add(bindingSelect, BorderLayout.CENTER);
            
        } catch (SQLException e) {
            e.printStackTrace();
        } catch (UnknownRoleException e) {
            e.printStackTrace();
        }

        // Options tab
        optionsPanel.setLayout(new GridLayout(4,1));
        relative = new JCheckBox("Relative vertical scale");
        hash = new JCheckBox("Show chromosome coordinates");
        common = new JCheckBox("Common vertical scale");
        seqletters = new JCheckBox("Show sequence");
        oldchipseq = new JCheckBox("Use old ChipSeq painter");
        hash.setSelected(true);
        optionsPanel.add(hash);
        optionsPanel.add(seqletters);
        optionsPanel.add(relative);
        optionsPanel.add(common);
        optionsPanel.add(oldchipseq);
        
        // Annotations tab
        JPanel lists = new JPanel();
        lists.setLayout(new GridLayout(3,2));
        JPanel boxes = new JPanel();
        boxes.setLayout(new GridLayout(3,2));
        genesmodel = new DefaultListModel();
        ncrnasmodel = new DefaultListModel();
        otherfeatsmodel = new DefaultListModel();
        genes = new JList(genesmodel);
        ncrnas = new JList(ncrnasmodel);
        otherfeats = new JList(otherfeatsmodel);
        genes.setVisibleRowCount(7);genes.setLayoutOrientation(JList.VERTICAL);
        otherfeats.setVisibleRowCount(7); otherfeats.setLayoutOrientation(JList.VERTICAL);
        ncrnas.setVisibleRowCount(7); ncrnas.setLayoutOrientation(JList.VERTICAL);
        
        polyA = new JCheckBox("PolyA sequences");
        gccontent = new JCheckBox("GC content");
        pyrpurcontent = new JCheckBox("Pyr/Pur content");
        cpg = new JCheckBox("CpG");
        regexmatcher = new JCheckBox("Regex Matcher");
        geneslabel = new JLabel("Genes");
        ncrnaslabel = new JLabel("ncRNAs");
        otherfeatslabel = new JLabel("Other annotations");

        lists.add(geneslabel);
        lists.add(new JScrollPane(genes));
        lists.add(ncrnaslabel);
        lists.add(new JScrollPane(ncrnas));
        lists.add(otherfeatslabel);
        lists.add(new JScrollPane(otherfeats));
        boxes.add(polyA);
        boxes.add(gccontent);
        boxes.add(pyrpurcontent);
        boxes.add(cpg);
        boxes.add(regexmatcher);

        annotationsPanel.setLayout(new BorderLayout());
        annotationsPanel.add(lists,BorderLayout.CENTER);
        annotationsPanel.add(boxes,BorderLayout.SOUTH);

        // file tracks tab
        filetracks = new FileBasedTracksPanel();        

        // use this to make the spacing be not-stupid
        JPanel dummy = new JPanel();        
        dummy.add(speciesLocationPanel);
        dummy.add(new JPanel());
        addTab("Species & Location",new JScrollPane(dummy));
        
        //dummy = new JPanel();  dummy.add(chipChipPanel); dummy.add(new JPanel());
        addTab("ChIP-Chip Data",chipChipPanel);

        addTab("Peaks",peakPanel);
        addTab("Expression", exprPanel);
        
        addTab("ChIP-Seq", chipSeqPanel);
        addTab("Paired ChIP-Seq", pairedChipSeqPanel);
        addTab("ChipSeq Analysis", chipSeqAnalysisPanel);


        dummy = new JPanel();  dummy.add(annotationsPanel); dummy.add(new JPanel());
        addTab("Annotations",new JScrollPane(dummy));

        dummy = new JPanel();  dummy.add(filetracks); dummy.add(new JPanel());
        addTab("File Tracks",new JScrollPane(dummy));
                
        dummy = new JPanel();  
        dummy.setLayout(new BorderLayout());
        dummy.add(motifScanPanel,BorderLayout.CENTER); 
        //        motifScanPanel.setPreferredSize(dummy.getPreferredSize());
        //        motifScanPanel.resetToPreferredSizes();
        //        addTab("Motifs",new JScrollPane(dummy));
        addTab("Motif Scans",dummy);

        dummy = new JPanel();
        dummy.setLayout(new BorderLayout());
        dummy.add(motifPanel,BorderLayout.CENTER);
        addTab("Motifs",dummy);
        
        dummy = new JPanel();  dummy.add(optionsPanel); dummy.add(new JPanel());
        addTab("Display Options",new JScrollPane(dummy));
    }

    public void init(WarpOptions opts) throws NotFoundException {
        handlingChange = true;
        init();
        createdFrom = opts;      
        
        if (opts.genome != null) {
            this.species.removeItemListener(this);
            this.genome.removeItemListener(this);
            this.species.setSelectedItem(opts.species);
            updateGenomeSelection();
            this.genome.setSelectedItem(opts.genome);            
            updateExptSelection();
            this.genome.addItemListener(this);
            this.species.addItemListener(this);
        } else if (opts.species != null) {
            this.species.setSelectedItem(opts.species);
            updateGenomeSelection();
        } else {
            this.species.setSelectedIndex(0);
            updateGenomeSelection();
        }
        handlingChange = false;        
        if (opts.gene != null &&
            !opts.gene.equals("")) {
            gene.setText(opts.gene);
        }
        relative.setSelected(opts.relative);
        hash.setSelected(opts.hash);
        gccontent.setSelected(opts.gccontent);        
        pyrpurcontent.setSelected(opts.pyrpurcontent);
        cpg.setSelected(opts.cpg);
        regexmatcher.setSelected(opts.regexmatcher);
        seqletters.setSelected(opts.seqletters);
        oldchipseq.setSelected(!opts.chipseqHistogramPainter);

        int[] selected = new int[opts.genes.size()];
        for (int i = 0; i < opts.genes.size(); i++) {
            selected[i] = genesmodel.indexOf(opts.genes.get(i));            
        }
        genes.setSelectedIndices(selected);
        selected = new int[opts.ncrnas.size()];
        for (int i = 0; i < opts.ncrnas.size(); i++) {
            selected[i] = ncrnasmodel.indexOf(opts.ncrnas.get(i));            
        }
        ncrnas.setSelectedIndices(selected);
        selected = new int[opts.otherannots.size()];
        for (int i = 0; i < opts.otherannots.size(); i++) {
            selected[i] = otherfeatsmodel.indexOf(opts.otherannots.get(i));
        }
        otherfeats.setSelectedIndices(selected);
        motifScanPanel.addToSelected(opts.motifscans);
        motifPanel.addToSelected(opts.motifs);
        
        try {
            Genome g = loadGenome();
            if (g != null) {
                ChipChipDataset ds = loadGenome().getChipChipDataset();
                for (int i = 0; i < opts.agilentdata.size(); i++) {
                    exptSelect.addToSelected(new ChipChipLocator(ds,
                                                                 opts.agilentdata.get(i).name,
                                                                 opts.agilentdata.get(i).version,
                                                                 opts.agilentdata.get(i).replicate));
                }
                for (int i = 0; i < opts.bayesresults.size(); i++) {
                    exptSelect.addToSelected(new BayesLocator(ds,
                                                              opts.bayesresults.get(i).name,
                                                              opts.bayesresults.get(i).version));
                }
                for (int i = 0; i < opts.msp.size(); i++) {
                    exptSelect.addToSelected(new MSPLocator(ds,
                                                            opts.msp.get(i).name,
                                                            opts.msp.get(i).version));
                }
            }
        } catch (NullPointerException ex) {
            /* this doesn't work if we can't get a genome object.  Just ignore
               the exception.  It only means that we can't fill in the
               selected experiments and teh user will have to do it again */
            ex.printStackTrace();
        }
        if (opts.position != null &&
            !opts.position.equals("")) {
            position.setText(opts.position);
        } else {
            if (opts.chrom != null &&
                !opts.chrom.equals("")) {
                position.setText(opts.chrom + ":" + opts.start + "-" + opts.stop);
            }
        }       
        chipSeqSelect.addToSelected(opts.chipseqExpts);
        pairedChipSeqSelect.addToSelected(opts.pairedChipseqExpts);
        bindingSelect.addToSelected(opts.bindingScans);
        exprSelect.addToSelected(opts.exprExperiments);
        filetracks.fill(opts.regionTracks);
    }

    public void setSpeciesGenomeDefaults() {
        String species = null, genome = null;
        try {
            ResourceBundle res = ResourceBundle.getBundle("edu.mit.csail.cgs.warpdrive.defaultSpeciesGenome");
            species = res.getString("species");
            genome = res.getString("genome");
        } catch (Exception e) {
            // don't do anything.  If it fails, then we just fall back to the defaults.
        }
        if (species == null || genome == null) {
            try {
                String homedir = System.getenv("HOME");
                String basename = "defaultSpeciesGenome";
                String fname = homedir + "/." + basename;
                File propfile = new File(fname);
                if (!(propfile.exists() && propfile.canRead())) {
                    homedir = System.getProperty("user.dir");
                    fname = homedir + "/" + basename;
                    propfile = new File(fname);
                }
                if (propfile.exists() && propfile.canRead()) {
                    InputStream is = new FileInputStream(propfile);
                    BufferedReader reader = new BufferedReader(new InputStreamReader(is));        
                    String line;
                    while ((line = reader.readLine()) != null) {
                        int p = line.indexOf('=');
                        String key = line.substring(0,p);
                        String value = line.substring(p+1);
                        if (key.equals("species")) {
                            species = value;
                        }
                        if (key.equals("genome")) {
                            genome = value;
                        }
                    }
                }
            } catch (Exception e) {
                // don't do anything.  If it fails, then we just fall back to the defaults.
            }

        }
        if (species == null || genome == null) {
            species = (String)this.species.getSelectedItem();
            genome = (String)this.genome.getSelectedItem();
        }
        this.species.setSelectedItem(species);
        this.genome.setSelectedItem(genome);
    }

    /* fills in and returns a WarpOptions object based on the current selections 
     */
    public WarpOptions parseOptions() {
        WarpOptions these = new WarpOptions();
        // parse the species and location tab
        these.species = species.getSelectedItem().toString();
        these.genome = genome.getSelectedItem().toString();
        these.position = position.getText();
        these.gene = gene.getText();

        // parse the options tab
        these.hash = hash.isSelected();
        these.relative = relative.isSelected();
        these.seqletters = seqletters.isSelected();
        these.chipseqHistogramPainter = !oldchipseq.isSelected();

        // parse the annotations tab
        Object[] selected = genes.getSelectedValues();
        for (int i = 0; i < selected.length; i++) {
            these.genes.add(selected[i].toString());
        }
        selected = ncrnas.getSelectedValues();
        for (int i = 0; i < selected.length; i++) {
            these.ncrnas.add(selected[i].toString());
        }
        selected = otherfeats.getSelectedValues();
        for (int i = 0; i < selected.length; i++) {
            these.otherannots.add(selected[i].toString());
        }
        these.motifscans.addAll(motifScanPanel.getSelected());
        these.motifs.addAll(motifPanel.getSelected());
        these.gccontent = gccontent.isSelected();
        these.pyrpurcontent = pyrpurcontent.isSelected();
        these.cpg = cpg.isSelected();
        these.regexmatcher = regexmatcher.isSelected();
        
        // parse the peaks tab
        Collection<BindingScan> scans = bindingSelect.getSelected();
        these.bindingScans.addAll(scans);
        
        // parse the expression tab
        Collection<Experiment> expts = exprSelect.getSelected();
        these.exprExperiments.addAll(expts);
        
        for(ChipSeqLocator loc : chipSeqSelect.getSelected()) { 
            these.chipseqExpts.add(loc);
        }
        for(ChipSeqLocator loc : pairedChipSeqSelect.getSelected()) { 
            these.pairedChipseqExpts.add(loc);
        }

        
        // parse the exptSelect panel selections.
        for(ExptLocator loc : exptSelect.getSelected()) { 
            
            if(loc instanceof ChipChipLocator) { 
                ChipChipLocator aloc = (ChipChipLocator)loc;
                these.agilentdata.add(aloc);
            }
            
            if(loc instanceof BayesLocator) { 
                BayesLocator bloc = (BayesLocator)loc;
                these.bayesresults.add(bloc);
            }
            
            if(loc instanceof MSPLocator) {
                MSPLocator mloc = (MSPLocator)loc;
                these.msp.add(mloc);
            }
        }
        filetracks.parse(these.regionTracks);
        these.regexes = regexes;
        return these;
    }
    
    public WarpOptions parseAndDiff() {
        WarpOptions these = parseOptions();
        // need to see if we have existing options and if they're compatible.
        // if they are, return the difference.  Otherwise, return the complete
        // options.
        if (createdFrom != null &&
            these.species.equals(createdFrom.species) &&
            these.genome.equals(createdFrom.genome)) {
            these.differenceOf(createdFrom);
        }
        return these;
    }
    
    /* updates the choice of experiments based on the
       currently selected genome and species */
    private void updateExptSelection() {
        Genome lg = loadGenome();
        Genome g = lg;
        
        System.err.println("UPDATING GENOME FOR EXPERIMENT SELECTION " + g);

        exptSelect.setGenome(lg);
        chipSeqSelect.setGenome(lg);
        pairedChipSeqSelect.setGenome(lg);
        motifPanel.setGenome(lg);
        motifScanPanel.setGenome(lg);
        if(lg != null) { 
            bindingSelect.setGenome(lg);            
        }
        // update the set of Gene annotations
        genesmodel.clear();
        otherfeatsmodel.clear();

        if(g != null) { 
            for(String type : gfLoader.getTypes(g)) {
                genesmodel.addElement(type);
            }
            for(String type : annotLoader.getTypes(g)) {
                otherfeatsmodel.addElement(type);
            }
            java.util.List<String> chroms = g.getChromList();
            if (chroms.size() == 0) {
                throw new RuntimeException("EMPTY CHROMOSOME LIST for " + g);
            }

            position.setText(chroms.get(0) + ":10000-20000");
        }
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
    
    public void actionPerformed (ActionEvent e) {

    }
}
