package edu.mit.csail.cgs.tools.chipchip;

import edu.mit.csail.cgs.utils.parsing.textfiles.*;
import edu.mit.csail.cgs.datasets.species.*;
import java.io.*;
import java.sql.*;
import java.util.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.Sequence;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.utils.parsing.alignment.PSL;
import edu.mit.csail.cgs.utils.parsing.alignment.PSLHit;
import edu.mit.csail.cgs.utils.parsing.textfiles.*;


/**
 * AddGAL (named for the old array design format) can load array 
 * designs in several formats (GAL, TDT, NDF).  You will need
 * a PSL file that maps the probes to the relevant genome.
 *
 *
 * For example:
 * AddGAL --species "$CE;Ce2" --design "2006-07-18" --designfile 2006-07-18_C_elegans_ChIP02.ndf --pslfile "probes.against_ce2.psl" 
 * 
 * If the design file does not contain the probe sequences, you can specify a FASTA file with --seqfile.
 *
 * You can also limit the set of genomic locations to which a probe is mapped with the
 *  --minmatch or --minscore options.  --minmatch is the number of basepairs that must match at a position
 *  and --minscore is the minimum score (computed as 2 * matches - mismatches).
 *
 * If an array design contains multiple slides, simply run this program with the same --design parameter once per slide.
 *
 * The schema keeps track of "galfile", which is parsed out of the --designfile name.  Since there are sometimes variants
 * on a name, ChipChipMetadataLoader.loadGALFile(String) does some parsing and conversion.  You may need to add
 * code there if you load an array design but none of its probes are found later when you try to load the data.
 *
 */
public class AddGAL implements ActionListener, ItemListener {
    
    private String species, genomeversion, designname;
    private ArrayList<String> designfile, pslfile, seqfile;
    private int speciesid, genomeid, designid;
    private ChipChipMetadataLoader chipchip;
    private java.sql.Connection cxn;
    private PreparedStatement insert, insertpos;
    private Pattern pospattern;
    private Genome genome;
    private Organism organism;
    private int posadded, minmatches, minscore;

    /* these are only used by the GUI code at the end of the file */
    private JComboBox speciesbox, genomebox, designbox;
    private JTextField minmatchesfield, minscorefield;        
    private JTextField designfield, pslfield, seqfield;
    private JButton addFileButton, removeFileButton, findDesignButton, findPSLButton, findSeqButton, okButton, cancelButton;
    private JList designlist, psllist, seqlist;
    private DefaultListModel designmodel, pslmodel, seqmodel;
    private AddGALTableModel tablemodel;
    private JTable table;
    private JFrame frame;

    public static void main(String args[]) {
        try {
            AddGAL adder = new AddGAL();
            if (args.length == 0) {
                adder.parseArgsGUI();
                return;
            } else {
                adder.parseArgs(args);
                adder.go();
                adder.cxn.commit();
            }
            System.err.println("Think I added " + adder.posadded + " positions");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public AddGAL() throws SQLException {
        designfile = new ArrayList<String>();
        pslfile = new ArrayList<String>();
        seqfile = new ArrayList<String>();
        pospattern = Pattern.compile("(.+)\\:(\\d+)\\-(\\d+)");
        posadded = 0;
        chipchip = new ChipChipMetadataLoader();            
        cxn = edu.mit.csail.cgs.utils.database.DatabaseFactory.getConnection("chipchip");
        cxn.setAutoCommit(false);
    }
    public void parseArgs(String args[]) throws NotFoundException, SQLException, UnknownRoleException {
        minmatches = -1;
        minscore = -1;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                String s[] = args[++i].split(";");
                species = s[0];
                if (s.length >= 2) {
                    genomeversion = s[1];
                }                
            }
            if (args[i].equals("--genomeversion")) {
                genomeversion  =args[++i];
            }
            if (args[i].equals("--design")) {
                designname = args[++i];
            }
            if (args[i].equals("--designfile")) {
                designfile.add(args[++i]);
            }
            if (args[i].equals("--pslfile")) {
                pslfile.add(args[++i]);
            }
            if (args[i].equals("--seqfile")) {
                seqfile.add(args[++i]);
            }
            if (args[i].equals("--minmatch") ||
                args[i].equals("--minmatches")) {
                minmatches = Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--minscore")) {
                minscore = Integer.parseInt(args[++i]);
            }
        }
        organism = new Organism(species);
        genome = organism.getGenome(genomeversion);
        speciesid = organism.getDBID();
        genomeid = genome.getDBID();
        designid = chipchip.loadArrayDesign(designname,genome).getDBID();        
    }    
    public void go() throws SQLException, IOException {
        insertpos = cxn.prepareStatement("insert into probelocation (id, chromosome, startpos, stoppos, strand,loccount, bitscore) values (?,?,?,?,?,?,?)");

        for (int i = 0; i < designfile.size(); i++) {
            String fname = designfile.get(i);
            String type = getType(fname);
            String name = fixName(fname);
            edu.mit.csail.cgs.datasets.chipchip.GALFile dbfile = 
                chipchip.loadGALFile(name);
            insert = cxn.prepareStatement("insert into probedesign(id, arraydesign, blockno, colno, rowno, galfile, probename, probeid, type, sequence) " + 
                                          "values("+Sequence.getInsertSQL(cxn,"probedesign_id")+"," + designid +",?,?,?," + dbfile.getDBID() + ",?,?,?,?)");        
            Map<String,ArrayList<Integer>> probes = new HashMap<String,ArrayList<Integer>>();
            Statement stmt = cxn.createStatement();
            ResultSet rs;
            RowsAndColumns file = null;
            RowsColumnsHandler handler;
            Map<String,String> probeseqs = null;
            if (i < seqfile.size()) {
                FASTAStream stream = new FASTAStream(new File(seqfile.get(i)));
                probeseqs = new HashMap<String,String>();
                while (stream.hasNext()) {
                    Pair<String,String> pair = stream.next();
                    probeseqs.put(pair.getFirst(), pair.getLast());
                }
                stream.close();
            }
            if (type.equals("GAL")) {
                edu.mit.csail.cgs.datasets.chipchip.AddGALHandler galhandler = new edu.mit.csail.cgs.datasets.chipchip.AddGALHandler(cxn,insert,probeseqs);
                file = new edu.mit.csail.cgs.utils.parsing.textfiles.GALFile(fname,galhandler);
                galhandler.setExistingKeys(galhandler.queryExistingKeys(designid,dbfile.getDBID()));
                handler = galhandler;
            } else if (type.equals("TDT")) {
                AddTDTHandler tdthandler = new AddTDTHandler(cxn,insert,probeseqs);
                file = new edu.mit.csail.cgs.utils.parsing.textfiles.TDTFile(fname,tdthandler);
                tdthandler.setExistingKeys(tdthandler.queryExistingKeys(designid,dbfile.getDBID()));
                handler = tdthandler;
            } else if (type.equals("NDF")) {
                AddNDFHandler ndfhandler = new AddNDFHandler(cxn,insert,probeseqs);
                file = new edu.mit.csail.cgs.utils.parsing.textfiles.NDFFile(fname,ndfhandler);
                ndfhandler.setExistingKeys(ndfhandler.queryExistingKeys(designid,dbfile.getDBID()));
                handler = ndfhandler;
            } else {
                throw new IllegalArgumentException("Didn't get back a known type " + type);
            }
            file.parse();           
            cxn.commit();
            probes.clear();
            rs = stmt.executeQuery("select probeid, id from probedesign where arraydesign = " + designid + 
                                   " and galfile = " + dbfile.getDBID());
            while (rs.next()) {
                if (!probes.containsKey(rs.getString(1))) {
                    probes.put(rs.getString(1),new ArrayList<Integer>());
                }
                probes.get(rs.getString(1)).add(rs.getInt(2));
            }
            rs.close(); stmt.close();            
            parseLocations(pslfile.get(i), probes);
            cxn.commit();
        }
    }
    public void parseLocations(String pslfname, Map<String,ArrayList<Integer>> probes) throws SQLException, IOException, FileNotFoundException  {
        Set<String> seen = new HashSet<String>();
        Map<String,Integer> chrommap = new HashMap<String,Integer>();
        Map<String,Integer> chromids = genome.getChromIDMap();
        BufferedReader reader = new BufferedReader(new FileReader(pslfname));
        PSL pslfile = new PSL(reader);
        String lastprobe = null;
        ArrayList<PSLHit> hits = new ArrayList<PSLHit>();
        while (pslfile.hasNext()) {
            PSLHit hit = pslfile.next();
            if (minmatches > -1) {
                if (hit.match < minmatches) {
                    continue;
                }
            }
            if (minscore > -1) {
                if (hit.match * 2 - hit.mismatch < minscore) {
                    continue;
                }
            }

            if (lastprobe != null &&
                !lastprobe.equals(hit.qname)) {
                if (!seen.contains(lastprobe)) {
                    addHits(hits,probes,chrommap,chromids);
                }
                hits.clear();
                seen.add(lastprobe);
            }
            hits.add(hit);
            lastprobe = hit.qname;
        }
        reader.close();
        addHits(hits,probes,chrommap,chromids);
    }
    public void addHits(ArrayList<PSLHit> hits, Map<String,ArrayList<Integer>> probes, Map<String,Integer> chrommap, Map<String,Integer> chromids) throws SQLException {
        for (PSLHit hit : hits) {
            if (!chrommap.containsKey(hit.tname)) {
                int id;
                String newchrom = hit.tname.replaceFirst("\\.fa.{0,3}$","").replaceFirst("^chr","");
                if (newchrom.equals("Mito")) {newchrom = "mt";}
                if (chromids.containsKey(newchrom)) {
                    id = chromids.get(newchrom);
                    chrommap.put(hit.tname,id);
                } else {
                    for (String s : chromids.keySet()) {
                        System.err.println("Knew about " + s);
                    }
                    throw new IllegalArgumentException("Can't figure out " + hit.tname + "," + newchrom);
                }
            }
            int chromid = chrommap.get(hit.tname);
            if (!probes.containsKey(hit.qname)) {continue;}
            for (int pid : probes.get(hit.qname)) {
                posadded++;
                insertpos.setInt(1,pid);
                insertpos.setInt(2,chromid);
                insertpos.setInt(3,hit.tstart);
                insertpos.setInt(4,hit.tend);
                insertpos.setString(5,Character.toString(hit.strand));
                insertpos.setInt(6,hits.size());
                insertpos.setDouble(7,2*hit.match - hit.mismatch);
                try {
                    insertpos.execute();
                } catch (SQLException e) {
                    if (e.getErrorCode() != 1) {
                        throw e;
                    }
                }
            }
            
        }
    }
    public static String getType(String fname) {
        if (fname.matches(".*gal$")) {
            return "GAL";
        } else if (fname.matches(".*tdt$")) {
            return "TDT";
        } else if (fname.matches(".*ndf$")) {
            return "NDF";
        } else {
            throw new IllegalArgumentException("Can't figure out file type of " + fname);
        }
    }
    public static String fixName (String fname) {
        fname = fname.replaceFirst("^.*/","");
        Pattern pattern = Pattern.compile("^(\\d+)_D.*");
        Matcher matcher = pattern.matcher(fname);
        if (matcher.matches()) {
            fname = matcher.group() + ".tdt";
            return fname;
        }
        return fname;
    }    

    /* graphical interface to adding an array design.
       Not only does this gather input, but it needs to arrange
       for this.go() and this.cxn.commit() to be called */
    public void parseArgsGUI() {
        frame = new JFrame();
        JPanel panel = new JPanel();
        Dimension buttonSize = new Dimension(30,20);
        GridBagLayout gridbag = new GridBagLayout();
        panel.setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();        
        constraints.weightx = 1.0;
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = GridBagConstraints.RELATIVE;
        JLabel label = new JLabel("Species");
        label.setHorizontalAlignment(SwingConstants.RIGHT);
        gridbag.setConstraints(label,constraints);
        panel.add(label);
        speciesbox = new JComboBox();
        speciesbox.addItemListener(this);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(speciesbox,constraints);
        panel.add(speciesbox);
        
        constraints.gridwidth = GridBagConstraints.RELATIVE;
        label = new JLabel("Genome");
        label.setHorizontalAlignment(SwingConstants.RIGHT);
        gridbag.setConstraints(label,constraints);
        panel.add(label);
        genomebox = new JComboBox();
        genomebox.addItemListener(this);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(genomebox,constraints);
        panel.add(genomebox);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        label = new JLabel("Array Design");
        label.setHorizontalAlignment(SwingConstants.RIGHT);
        gridbag.setConstraints(label,constraints);
        panel.add(label);
        designbox = new JComboBox();
        designbox.addItemListener(this);
        designbox.setEditable(true);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(designbox,constraints);
        panel.add(designbox);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        label = new JLabel("Minimum base matches");
        label.setHorizontalAlignment(SwingConstants.RIGHT);
        gridbag.setConstraints(label,constraints);
        panel.add(label);
        minmatchesfield = new JTextField("45");
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(minmatchesfield,constraints);
        panel.add(minmatchesfield);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        label = new JLabel("Minimum score");
        label.setHorizontalAlignment(SwingConstants.RIGHT);
        gridbag.setConstraints(label,constraints);
        panel.add(label);
        minscorefield = new JTextField("45");
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(minscorefield,constraints);
        panel.add(minscorefield);

        tablemodel = new AddGALTableModel();
        table = new JTable(tablemodel);
        table.getTableHeader().setResizingAllowed(true);
        table.getTableHeader().getColumnModel().getColumn(0).setHeaderValue("Design File");
        table.getTableHeader().getColumnModel().getColumn(1).setHeaderValue("PSL File");
        table.getTableHeader().getColumnModel().getColumn(2).setHeaderValue("Sequence File");
        gridbag.setConstraints(table,constraints);
        //        JScrollPane sp = new JScrollPane();
        //        sp.add(table);
        //        panel.add(sp);
        panel.add(table);

        addFileButton = new JButton("Add Files");
        removeFileButton = new JButton("Remove Files");
        addFileButton.addActionListener(this);
        removeFileButton.addActionListener(this);
        constraints.gridwidth = GridBagConstraints.RELATIVE;
        gridbag.setConstraints(addFileButton,constraints);
        panel.add(addFileButton);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(removeFileButton,constraints);
        panel.add(removeFileButton);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        findDesignButton = new JButton("Find Design");
        designfield = new JTextField();
        gridbag.setConstraints(designfield,constraints);
        panel.add(designfield);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(findDesignButton,constraints);
        panel.add(findDesignButton);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        findPSLButton = new JButton("Find PSL");
        pslfield = new JTextField();
        gridbag.setConstraints(pslfield,constraints);
        panel.add(pslfield);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(findPSLButton,constraints);
        panel.add(findPSLButton);

        constraints.gridwidth = GridBagConstraints.RELATIVE;
        findSeqButton = new JButton("Find Seq");
        seqfield = new JTextField();
        gridbag.setConstraints(seqfield,constraints);
        panel.add(seqfield);
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(findSeqButton,constraints);
        panel.add(findSeqButton);

        findDesignButton.addActionListener(this);
        findPSLButton.addActionListener(this);
        findSeqButton.addActionListener(this);
        okButton = new JButton("OK");
        cancelButton = new JButton("Cancel");
        okButton.addActionListener(this);
        cancelButton.addActionListener(this);        
        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1,2));
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);

        findDesignButton.setMaximumSize(buttonSize);
        findPSLButton.setMaximumSize(buttonSize);
        findSeqButton.setMaximumSize(buttonSize);
        okButton.setMaximumSize(buttonSize);
        cancelButton.setMaximumSize(buttonSize);

        setSpecies();
        setGenomes();
        JPanel outer = new JPanel();
        outer.setLayout(new BorderLayout());
        outer.add(panel,BorderLayout.CENTER);
        outer.add(buttonPanel,BorderLayout.SOUTH);
        frame.setContentPane(outer);
        frame.setSize(500,650);
        frame.setLocation(50,50);
        frame.setVisible(true);
    }

    public void setSpecies() {
        Collection<String> names = Organism.getOrganismNames();
        speciesbox.removeAllItems();
        for (String s : names) {
            speciesbox.addItem(s);
        }
        speciesbox.setSelectedIndex(0);

    }
    public void setGenomes() {
        try {
            String speciesname = (String)speciesbox.getSelectedItem();
            Organism o = new Organism(speciesname);
            Collection<String> names = o.getGenomeNames();
            genomebox.removeAllItems();
            for (String s : names) {
                genomebox.addItem(s);
            }
            designbox.removeAllItems();
            for (ArrayDesign d : chipchip.loadAllArrayDesigns(o)) {
                designbox.addItem(d.getName());
            }
        } catch (NotFoundException e) {
            System.err.println("Couldn't find species " + speciesbox.getSelectedItem());
            e.printStackTrace();
        } catch (SQLException e) {
            System.err.println("Couldn't find species or genomes" + speciesbox.getSelectedItem());
            e.printStackTrace();
        }

    }

    public void actionPerformed (ActionEvent e) {
        Object source = e.getSource();

        if (source == findDesignButton) {
            JFileChooser chooser = new JFileChooser(new File(System.getProperty("user.dir")));
            chooser.setDialogTitle("Select Array Design Files");
            chooser.setDialogType(JFileChooser.SAVE_DIALOG);
            chooser.setMultiSelectionEnabled(false);
            if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
                designfield.setText(chooser.getSelectedFile().getPath());
            }
        }

        if (source == findPSLButton) {
            JFileChooser chooser = new JFileChooser(new File(System.getProperty("user.dir")));
            chooser.setDialogTitle("Select PSL Files");
            chooser.setDialogType(JFileChooser.SAVE_DIALOG);
            chooser.setMultiSelectionEnabled(false);
            if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
                pslfield.setText(chooser.getSelectedFile().getPath());
            }
        }

        if (source == findSeqButton) {
            JFileChooser chooser = new JFileChooser(new File(System.getProperty("user.dir")));
            chooser.setDialogTitle("Select Probe Sequence Files");
            chooser.setDialogType(JFileChooser.SAVE_DIALOG);
            chooser.setMultiSelectionEnabled(false);
            if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
                seqfield.setText(chooser.getSelectedFile().getPath());
            }
        }

        if (source == addFileButton) {
            tablemodel.addFiles(designfield.getText(),
                                pslfield.getText(),
                                seqfield.getText());
        }

        if (source == removeFileButton) {
            int row = table.getSelectedRow();
            if (row >= 0) {
                tablemodel.removeFile(row);
            }
        }
        
        if (source == okButton) {
            try {
                organism = new Organism((String)speciesbox.getSelectedItem());
                genome = organism.getGenome((String)genomebox.getSelectedItem());
                speciesid = organism.getDBID();
                genomeid = genome.getDBID();
                designid = chipchip.loadArrayDesign((String)designbox.getSelectedItem(),genome).getDBID();        
                for (int row = 0; row < tablemodel.getRowCount(); row++) {
                    designfile.add((String)tablemodel.getValueAt(row,0));
                    pslfile.add((String)tablemodel.getValueAt(row,1));
                    seqfile.add((String)tablemodel.getValueAt(row,2));                    
                }
                go();
                cxn.commit();
            } catch (NotFoundException ex) {
                ex.printStackTrace();
            } catch (SQLException  ex) {
                ex.printStackTrace();
            } catch (IOException  ex) {
                ex.printStackTrace();
            }


        }
        if (source == cancelButton) {
            frame.dispose();
        }



    }
    public void itemStateChanged(ItemEvent e) {
        Object source = e.getItemSelectable();
        if (source == speciesbox) {
            setGenomes();
        }        
    }
}

class AddGALTableModel extends DefaultTableModel {
    private ArrayList<String> design, psl, seq;

    public AddGALTableModel() {
        design = new ArrayList<String>();
        psl = new ArrayList<String>();
        seq = new ArrayList<String>();
    }
    public void addFiles(String d, String p, String s) {
        design.add(d);
        psl.add(p);
        seq.add(s);
        fireTableRowsInserted(design.size(),
                              design.size());
    }
    public void removeFile(int i) {
        if (i >= design.size()) {
            return;
        }
        design.remove(i);
        psl.remove(i);
        seq.remove(i);
        fireTableRowsDeleted(i,i);
    }
    public Class getColumnClass(int c) {
        try {
            return Class.forName("java.lang.String");
        } catch (ClassNotFoundException e) {
            throw new RuntimeException("Couldn't get class for java.lang.String",e);
        }
    }
    public int getColumnCount() {
        return 3;
    }
    public String getColumnName(int col) {
        if (col == 0) {
            return "Design File";
        } else if (col == 1) {
            return "PSL File";
        } else if (col == 2) {
            return "FASTA probe seqs";
        } else {
            return null;
        }

    }
    public int getRowCount() {
        int act = design == null ? 0 : design.size();
        if (act < 5) {
            return 5;
        } else {
            return act;
        }
    }

    public Object getValueAt(int row, int col) {
        if (row >= design.size()) {
            return null;
        }
        if (col == 0) {
            return design.get(row);
        } else if (col == 1) {
            return psl.get(row);
        } else if (col == 2) {
            return seq.get(row);
        } else {
            throw new ArrayIndexOutOfBoundsException ("No such column " + col);
        }

    }
    public boolean isCellEditable(int row, int col) {return false;}
}

