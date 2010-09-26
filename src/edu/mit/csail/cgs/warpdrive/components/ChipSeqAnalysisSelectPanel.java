package edu.mit.csail.cgs.warpdrive.components;

import java.util.*;
import java.util.regex.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.sql.*;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.viz.components.GenericSelectPanel;

public class ChipSeqAnalysisSelectPanel extends GenericSelectPanel<ChipSeqAnalysis> {

    TreeSet<ChipSeqAnalysis> analyses;
    private JTextField regex;
    private JCheckBox active;
    private ChipSeqAnalysisTableModel selectedModel, filteredModel;

    public ChipSeqAnalysisSelectPanel(Genome g) { 
        super(g);
        analyses = new TreeSet<ChipSeqAnalysis>();
        selectedModel = new ChipSeqAnalysisTableModel();
        filteredModel = new ChipSeqAnalysisTableModel();
        init(filteredModel,selectedModel);
    }
    public ChipSeqAnalysisSelectPanel() { 
        super();
        analyses = new TreeSet<ChipSeqAnalysis>();
        selectedModel = new ChipSeqAnalysisTableModel();
        filteredModel = new ChipSeqAnalysisTableModel();
        init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());
        inputPanel.setLayout(new BorderLayout());
        regex = new JTextField();
        active = new JCheckBox("Active Analyses Only?", true);
        inputPanel.add(active, BorderLayout.EAST);
        inputPanel.add(new JLabel("pattern to filter analyses"), BorderLayout.WEST);
        inputPanel.add(regex, BorderLayout.CENTER);        
        return inputPanel;
    }
    
    public void retrieveData() {
        analyses.clear();
        try {
            synchronized(analyses) {
                Collection<ChipSeqAnalysis> all = ChipSeqAnalysis.getAll(active.isSelected() ? true : null);
                for(ChipSeqAnalysis a :all) { 
                    analyses.add(a);
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException(e.toString(), e);
        }
    }
    public void updateComponents() {
        selectedModel.clear();
        filteredModel.clear();
        synchronized(analyses) {
            for (ChipSeqAnalysis a : analyses) {
                filteredModel.addObject(a);
            }
        }
    }
    public void filter() {
        String reg = regex.getText().trim();
        Pattern patt = null;
        if(reg != null && reg.length() > 0) {
            patt = Pattern.compile(reg);
        }
        synchronized(analyses) {
            analyses.clear();
            try {
                Collection<ChipSeqAnalysis> all = ChipSeqAnalysis.getAll(active.isSelected() ? true : null);
                for (ChipSeqAnalysis a : all) {
                    if (patt == null || patt.matcher(a.toString()).find()) {
                        Set<ChipSeqAlignment> fg = a.getForeground();
                        Iterator<ChipSeqAlignment> i = fg.iterator();                        
                        if (i.hasNext()) {
                            ChipSeqAlignment align = i.next();
                            if (!align.getGenome().equals(getGenome())) {
                                continue;
                            }                            
                            analyses.add(a);
                        }
                    }
                }
            } catch (SQLException e) {
                throw new RuntimeException(e.toString(), e);
            }
            filteredModel.clear();
            for (ChipSeqAnalysis a : analyses) {
                filteredModel.addObject(a);
            }
        }
    }
    public void close() {
        super.close();
    }
}




