package edu.mit.csail.cgs.tools.chipseq;

import java.util.Collection;
import java.util.Collections;
import java.util.ArrayList;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;


public class AnalysisInputs {
    
    public static void main(String[] args) throws Exception {
        Collection<ChipSeqAnalysis> analyses = Args.parseChipSeqAnalyses(args,"chipseqanalysis");
        for (ChipSeqAnalysis analysis : analyses) {
            for (ChipSeqAlignment align : analysis.getForeground()) {
                ChipSeqExpt expt = align.getExpt();
                System.out.println(analysis.toString() + "\t" + "FOREGROUND" + "\t" + 
                                   expt.getName() + ";" + expt.getReplicate() + ";" + align.getName());
            }
            for (ChipSeqAlignment align : analysis.getBackground()) {
                ChipSeqExpt expt = align.getExpt();
                System.out.println(analysis.toString() + "\t" + "BACKGROUND" + "\t" + 
                                   expt.getName() + ";" + expt.getReplicate() + ";" + align.getName());
            }

        }

    }


}