package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.sql.*;
import java.io.*;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;

/** Generates a .png file with a set of motif logos in it.
 * Needs --species and --out (output file name) on the command line.
 * Reads motifs on STDIN.  Fields are tab or semicolon separated and 
 * are name, version, and an optional additional label
 *
 */

public class DrawMotifs {

    public static void main(String args[]) throws Exception {
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
        String outfilename = Args.parseString(args,"out","motifs.png");
        Genome genome = Args.parseGenome(args).cdr();
        int species = genome.getSpeciesDBID();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        WeightMatrixLoader loader = new WeightMatrixLoader();
        List<WeightMatrix> matrices = new ArrayList<WeightMatrix>();
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("[\\t;]");
            for (WeightMatrix m : loader.query(pieces[0], pieces[1], null)) {
                if (pieces.length >= 3) {
                    m.version = m.version + " (" + pieces[2] + ")";
                }
                matrices.add(m);
            }
        }
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }
        ClusterMotifs.drawCluster(matrices,
                                  outfilename,
                                  800,
                                  200,
                                  2);



    }

}