package edu.mit.csail.cgs.tools.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;

/**
 * returns the ids of the alignments specified on the command line.
 * All fields are optional.  factor, cells, and condition are ANDed together.  align 
 * parameters are ORed together and then applied to the factor/cells/condition.  The
 * fields include --align name;replicate;version
 *                --cells
 *                --factor
 *                --condition
 *
 * GetAlignments --species "$SC;SGDv1" --align "name;replicate;alignment version" --factor "Gcn4" --condition "YPD"
 */

public class GetAlignments {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipseq");
        cxn.setAutoCommit(false);
        Genome genome = Args.parseGenome(args).cdr();
        Collection<String> alignnames = Args.parseStrings(args,"align");
        String factorstring = Args.parseString(args,"factor",null);
        String cellsstring = Args.parseString(args,"cells",null);
        String conditionstring = Args.parseString(args,"condition",null);

        ChipSeqLoader loader = new ChipSeqLoader();
        MetadataLoader core = new MetadataLoader();

        Integer factor = null, cells = null, condition = null;
        if (factorstring != null) {
            factor = core.getFactor(factorstring).getDBID();
        }
        if (cellsstring != null) {
            cells = core.getCells(cellsstring).getDBID();
        }
        if (conditionstring != null) {
            condition = core.getCondition(conditionstring).getDBID();
        }

        for (String an : alignnames) {
            String pieces[] = an.split(";");
            for (ChipSeqAlignment a : loader.loadAlignments(pieces[0],
                                                            (pieces.length > 1 && pieces[1] != "") ? pieces[1] : null,
                                                            (pieces.length > 2 && pieces[2] != "") ? pieces[2] : null,
                                                            factor,cells,condition,genome)) {
                System.out.println(a.getDBID());
            }
        }
    }
}