package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import java.text.DecimalFormat;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.io.parsing.FASTAStream;
import edu.mit.csail.cgs.utils.*;
//import edu.mit.csail.cgs.utils.probability.Binomial;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.tools.motifs.*;
import edu.mit.csail.cgs.tools.utils.Args;

public class SVMFiles extends CombinatorialEnrichment {
    double trainfrac = .3;
    private String prefix;
    public SVMFiles() {
        super();
    }
    public void parseArgs(String args[]) throws Exception {
        super.parseArgs(args);
        trainfrac = Args.parseDouble(args,"trainfrac",trainfrac);
        prefix = Args.parseString(args,"prefix","svm");
    }
    private void saveLine(PrintWriter file, double val, WMHit[] hits) throws IOException {
        StringBuffer line = new StringBuffer(Double.toString(val));
        for (int i = 0; i < hits.length; i++) {
            line.append(String.format(" %d:%.4f",i, hits[i] == null ? 0 : hits[i].getScore()));
        }
        file.println(line.toString());
    }
    public void saveFiles() throws IOException {
        PrintWriter training = new PrintWriter(prefix + "_training.txt");
        PrintWriter test = new PrintWriter(prefix + "_test.txt");
        PrintWriter trainingregions = new PrintWriter(prefix + "_training.regions");
        PrintWriter testregions = new PrintWriter(prefix + "_test.regions");
        for (String s : fghits.keySet()) {
            if (Math.random() <= trainfrac) {
                saveLine(training, 1, fghits.get(s));
                trainingregions.println(s);
            } else {
                saveLine(test, 1,fghits.get(s));
                testregions.println(s);
            }
        }
        for (String s : bghits.keySet()) {
            if (Math.random() <= trainfrac) {
                saveLine(training, -1, bghits.get(s));
            } else {
                saveLine(test,-1, bghits.get(s));
                testregions.println(s);
            }
        }
        training.close();
        test.close();
        trainingregions.close();
        testregions.close();

    }

    public static void main(String args[]) throws Exception {
        SVMFiles ce = new SVMFiles();
        ce.parseArgs(args);
        System.err.println("Masking and saving");
        ce.maskSequence();
        ce.saveSequences();
        System.err.println("Doing weight matrix scanning");
        ce.doScans();
        System.err.println("Saving files");
        ce.saveFiles();
    }
        

}