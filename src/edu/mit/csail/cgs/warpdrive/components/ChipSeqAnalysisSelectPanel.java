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
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.viz.components.GenericSelectPanel;

public class ChipSeqAnalysisSelectPanel extends GenericSelectPanel<ChipSeqAnalysis> {

    TreeSet<ChipSeqAnalysis> analyses;
    private JTextField regex;
    private ChipSeqAnalysisTableModel selectedModel, filteredModel;

    public ChipSeqAnalysisSelectPanel() { 
        analyses = new TreeSet<ChipSeqAnalysis>();
        selectedModel = new ChipSeqAnalysisTableModel();
        filteredModel = new ChipSeqAnalysisTableModel();
        init(filteredModel,selectedModel);
    }
    public JPanel getInputsPanel() {
        JPanel inputPanel = new JPanel(); inputPanel.setLayout(new BorderLayout());
        inputPanel.setLayout(new BorderLayout());
        regex = new JTextField();
        inputPanel.add(new JLabel("pattern to filter analyses"), BorderLayout.WEST);
        inputPanel.add(regex, BorderLayout.CENTER);        
        return inputPanel;
    }
    
    public void retrieveData() {
        analyses.clear();
        try {
            synchronized(analyses) {
                Collection<ChipSeqAnalysis> all = ChipSeqAnalysis.getAll();
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
                for (ChipSeqAnalysis a : ChipSeqAnalysis.getAll()) {
                    Set<ChipSeqAlignment> fg = a.getForeground();
                    Iterator<ChipSeqAlignment> i = fg.iterator();
                    if (!i.next().getGenome().equals(getGenome())) {
                        continue;
                    }
                    if (patt == null || patt.matcher(a.toString()).find()) {
                        analyses.add(a);
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




