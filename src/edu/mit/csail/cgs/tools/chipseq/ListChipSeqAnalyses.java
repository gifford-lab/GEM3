package edu.mit.csail.cgs.tools.chipseq;

import java.util.Collection;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;

public class ListChipSeqAnalyses {

    public static void main(String[] args) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        Collection<ChipSeqAnalysis> all = ChipSeqAnalysis.getAll();
        for (ChipSeqAnalysis a : all) {
            int count = a.countResults(genome);
            if (count == 0) {
                continue;
            }
            System.err.println(String.format("%s\t%s\t%s\t%d\t%s\t%s",
                                             a.getName(),
                                             a.getVersion(),
                                             a.getProgramName(),
                                             count,
                                             a.getForeground().toString(),
                                             a.getBackground().toString()));
        }
    }

}