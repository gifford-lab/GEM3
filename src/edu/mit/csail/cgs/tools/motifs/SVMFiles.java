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
    private double motifWidth = 12; // width of the motif on which regions are centered
    private int kmerLength, maxOffset;
    private List<KmerFeature> features;
    public SVMFiles() {
        super();
    }
    public void parseArgs(String args[]) throws Exception {
        super.parseArgs(args);
        trainfrac = Args.parseDouble(args,"trainfrac",trainfrac);
        prefix = Args.parseString(args,"prefix","svm");
        kmerLength = Args.parseInt(args,"k",4);
        maxOffset = Args.parseInt(args,"maxoffset",4);

        features = new ArrayList<KmerFeature>();
        for (int i = 0; i < Math.pow(4,k); i++) {
            KmerFeature f = new KmerFeature();
            f.kmer = 
            
            f.offset = 0;
        }
        for (int pos = 0; pos < k; pos++) {
            int chunksize = Math.pow(4, (k - pos));

        }

        for (int o = 0; o < maxOffset; o++) {

        }

    }
    private void saveLine(PrintWriter file, double val, WMHit[] hits, double kmerFeatures[]) throws IOException {
        StringBuffer line = new StringBuffer(Double.toString(val));
        for (int i = 0; i < hits.length; i++) {
            line.append(String.format(" %d:%.4f",i+1, hits[i] == null ? 0 : hits[i].getScore()));
        }
        for (int i = 0; i < kmerFeatures.length; i++) {
            line.append(String.format(" %d:%.4f",i+hits.length+1, kmerFeatures[i]));
        }
        file.println(line.toString());
    }
    private double[] kmerFeatures(char[] sequence) {
        double output[] = new double[features.size()];
        for (int i = 0; i < features.size(); i++) {
            boolean match = true;
            KmerFeature feature = features.get(i);
            for (int j = 0; j < feature.kmer.length; j++) {
                if (sequence[feature.offset + j] != feature.kmer[j]) {
                    match = false;
                    break;
                }
            }
            output[i] = match ? 1 : 0;
        }
        return output;
    }
    public void saveFiles() throws IOException {
        PrintWriter training = new PrintWriter(prefix + "_training.txt");
        PrintWriter test = new PrintWriter(prefix + "_test.txt");
        PrintWriter trainingregions = new PrintWriter(prefix + "_training.regions");
        PrintWriter testregions = new PrintWriter(prefix + "_test.regions");
        PrintWriter featurenames = new PrintWriter(prefix + "_featurenames.txt");
        for (String s : fghits.keySet()) {
            if (Math.random() <= trainfrac) {
                saveLine(training, 1, fghits.get(s), kmerFeatures(foreground.get(s)));
                trainingregions.println(s);
            } else {
                saveLine(test, 1,fghits.get(s), kmerFeatures(foreground.get(s)));
                testregions.println(s);
            }
        }
        for (String s : bghits.keySet()) {
            if (Math.random() <= trainfrac) {
                saveLine(training, -1, bghits.get(s), kmerFeatures(background.get(s)));
            } else {
                saveLine(test,-1, bghits.get(s), kmerFeatures(background.get(s)));
                testregions.println(s);
            }
        }
        for (WeightMatrix m : matrices) {
            featurenames.println(m.toString());
        }

        training.close();
        test.close();
        trainingregions.close();
        testregions.close();
        featurenames.close();

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
class KmerFeature {
    public int offset;
    public char[] kmer;

}