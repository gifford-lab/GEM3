package edu.mit.csail.cgs.tools.chipseq;

import java.io.*;
import java.sql.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.alignments.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;

/** creates an experiment (if necessary) and alignment (if necessary)
 * in the database and prints the DBID on stdout.  Use this with
 * the readdb importer since we're still using the old oracle
 * stuff for metadata
 *
 * Usage:
 * CreateAlignment --species "$SC;SGDv1" --align "name;replicate;alignment version" --factor "Gcn4" --cells "FY4" --condition "YPD" --paramsfile params.txt --readlength 36
 */
public class CreateAlignment {

    public static void main(String args[]) throws SQLException, IOException, NotFoundException {

        java.sql.Connection cxn = DatabaseFactory.getConnection("chipseq");
        cxn.setAutoCommit(false);
        Genome genome = Args.parseGenome(args).cdr();
        String alignname = Args.parseString(args,"align",null);
        String alignpieces[] = alignname.split(";");
        String factorstring = Args.parseString(args,"factor",null);
        String cellsstring = Args.parseString(args,"cells",null);
        String conditionstring = Args.parseString(args,"condition",null);
        String paramsfname = Args.parseString(args,"paramsfile",null);
        int readlength = Args.parseInteger(args,"readlength",36);
        
        ChipSeqExpt expt = null;
        ChipSeqAlignment alignment = null;
        ChipSeqLoader loader = new ChipSeqLoader();
        MetadataLoader core = new MetadataLoader();
        try {
            expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
        } catch (NotFoundException e) {
            System.err.println("Creating experiment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
            PreparedStatement insert = ChipSeqExpt.createInsert(cxn);
            insert.setString(1, alignpieces[0]);
            insert.setString(2, alignpieces[1]);
            insert.setInt(3, genome.getSpeciesDBID());
            insert.setInt(4, readlength);
            insert.setInt(5, core.getCells(cellsstring).getDBID());
            insert.setInt(6, core.getCondition(conditionstring).getDBID());
            insert.setInt(7, core.getFactor(factorstring).getDBID());
            insert.execute();
            try {
                expt = loader.loadExperiment(alignpieces[0], alignpieces[1]);
            } catch (NotFoundException e2) {
                /* failed again means the insert failed.  you lose */
                cxn.rollback();
                throw new DatabaseException("Couldn't create " + alignpieces[0] + "," + alignpieces[1]);
            }
        }
        alignment = loader.loadAlignment(expt, alignpieces[2], genome);

        if (alignment == null) {
            try {
                PreparedStatement insert = ChipSeqAlignment.createInsertStatement(cxn);
                System.err.println("Creating alignment " + alignpieces[0] + ";" + alignpieces[1] + ";" + alignpieces[2]);
                System.err.println("Inserting alignment for experiment " + expt.getDBID());
                insert.setInt(1, expt.getDBID());
                insert.setString(2, alignpieces[2]);
                insert.setInt(3, genome.getDBID());
                insert.executeQuery();
                alignment = loader.loadAlignment(expt, alignpieces[2], genome);
                cxn.commit();
                File f = null;
                if (paramsfname != null) {
                    f = new File(paramsfname);
                }
                if (f != null && f.exists()) {
                    System.err.println("Reading alignment parameters from " + f);
                    loader.addAlignmentParameters(alignment, f);

                }
			} catch (IOException e) {
                cxn.rollback();
                System.err.println("Couldn't add alignment parameters");
                e.printStackTrace();
            }

        }
        loader.close();
        if (alignment == null) {
            cxn.rollback();
            throw new DatabaseException("Couldn't create alignment " + alignpieces[2] + " for " + alignpieces[0]);
        }
        System.out.println(alignment.getDBID());
        cxn.commit();
        
    }


}