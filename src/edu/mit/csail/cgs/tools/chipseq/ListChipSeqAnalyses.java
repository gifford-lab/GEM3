package edu.mit.csail.cgs.tools.chipseq;

import java.util.Collection;
import java.util.Collections;
import java.util.ArrayList;
import java.sql.SQLException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * output fields are
 * 0) name
 * 1) version
 * 2) program name
 * 3) event count
 * 4) fg expts
 * 5) bg expts
 * 6) dbid
 *
 */

public class ListChipSeqAnalyses {

    public static void main(String[] args) throws Exception {
        Genome genome = Args.parseGenome(args).cdr();
        boolean html = Args.parseFlags(args).contains("html");
        ArrayList<ChipSeqAnalysis> all = new ArrayList<ChipSeqAnalysis>();
        all.addAll(ChipSeqAnalysis.getAll());
        Collections.sort(all);
        if (html) {
            System.out.println("<tr><th>Name</th><th>Version</th><th>Analysis Program</th><th>Number of Binding Events</th><th>IP Experiments</th><th>Control Experimentss</th><th>database id</th></tr>");
        } 


        for (ChipSeqAnalysis a : all) {
            int count = a.countResults(genome);
            if (count == 0) {
                continue;
            }
            if (html) {
                System.out.println(String.format("<tr><td>%s</td><td>%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%d</td></tr>",
                                                 a.getName(),
                                                 a.getVersion(),
                                                 a.getProgramName(),
                                                 count,
                                                 a.getForeground().toString(),
                                                 a.getBackground().toString(),
                                                 a.getDBID()));
            } else {
                System.out.println(String.format("%s\t%s\t%s\t%d\t%s\t%s\t%d",
                                                 a.getName(),
                                                 a.getVersion(),
                                                 a.getProgramName(),
                                                 count,
                                                 a.getForeground().toString(),
                                                 a.getBackground().toString(),
                                                 a.getDBID()));
            }
        }
    }

}