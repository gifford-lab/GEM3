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

/** Generate input files for libsvm or liblinear classification.
 * Uses both motifs (see parameters for CompareEnrichment) and short kmers. 
 * The short kmers assume the regions are centered around a motif of size --motifwidth
 * 
 * Specify the kmers with
 * --k 4
 * --maxoffset 3 
 *
 * java SVMFiles --species "$MM;mm9" --first foo.txt --second bar.txt --k 4 --maxoffset 2 --motifwidth 12 --expand 6
 *
 * Where foo.fasta and bar.txt are files of genomic points (single base positions).  Or
 * you could use width 0 where the files are genomic regions of width 12.
 *
 * TODO:
 * - right now the kmer stuff doesn't look for reverse complements
 *
 */

public class SVMFiles extends CombinatorialEnrichment {
    double trainfrac = .3;
    private String prefix;
    private int kmerLength, maxOffset, motifWidth, basicKmers;
    private List<KmerFeature> features;
    private boolean stranded;
    public SVMFiles() {
        super();
    }
    public void parseArgs(String args[]) throws Exception {
        super.parseArgs(args);
        trainfrac = Args.parseDouble(args,"trainfrac",trainfrac);
        prefix = Args.parseString(args,"prefix","svm");
        kmerLength = Args.parseInteger(args,"k",4);
        maxOffset = Args.parseInteger(args,"maxoffset",4);
        motifWidth = Args.parseInteger(args,"motifwidth",12);
        stranded = Args.parseFlags(args).contains("stranded");

        features = new ArrayList<KmerFeature>();
        basicKmers = (int)(Math.pow(4,kmerLength));
        List<KmerFeature> basicFeatures = new ArrayList<KmerFeature>();
        for (int i = 0; i < basicKmers; i++) {
            KmerFeature f = new KmerFeature();
            f.kmer = DiscriminativeKmers.longToChars(i,kmerLength);
            f.offset = 0;
            basicFeatures.add(f);
        }
        System.err.println("There are " + basicKmers + " basic kmers");
        for (int o = 1; o <= maxOffset; o++) {
            for (int i = 0; i < basicKmers; i++) {
                KmerFeature basic = basicFeatures.get(i);                
                KmerFeature f = new KmerFeature();
                f.kmer = basic.kmer;
                f.offset = o;
                features.add(f);
            }
        }
        if (stranded) {
            for (int o = 1; o <= maxOffset; o++) {
                for (int i = 0; i < basicKmers; i++) {
                    KmerFeature basic = basicFeatures.get(i);                
                    KmerFeature f = new KmerFeature();
                    f.kmer = basic.kmer;
                    f.offset = -1 * o;
                features.add(f);
                }
            }
        }
        System.err.println("There are " + features.size() + " kmer features");
    }
    private void saveLine(PrintWriter file, double val, WMHit[] hits, double kmerFeatures[]) throws IOException {
        StringBuffer line = new StringBuffer(Double.toString(val));
        for (int i = 0; i < hits.length; i++) {
            line.append(String.format(" %d:%.4f",i+1, hits[i] == null ? 0 : hits[i].getScore()));
        }
        for (int i = 0; i < kmerFeatures.length; i++) {
            if (kmerFeatures[i] > 0) {
                line.append(String.format(" %d:%.4f",i+hits.length+1, kmerFeatures[i]));
            }
        }
        file.println(line.toString());
    }
    private double[] kmerFeatures(char[] sequence) {
        int sideSize = (sequence.length - motifWidth) / 2;
        double output[] = new double[features.size()];
        char[] tmp = new char[kmerLength];
        int offset = sideSize + motifWidth;
        for (int o = 1; o <= maxOffset; o++) {
            for (int j = 0; j < kmerLength; j++) {
                tmp[j] = sequence[offset+o+j];
            }
            long bk = DiscriminativeKmers.charsToLong(tmp);
            long ak = basicKmers * (o -1) + bk;
            output[(int)ak] = 1;
        }
        offset = sideSize - maxOffset - kmerLength;
        for (int o = 1; o <= maxOffset; o++) {
            for (int j = 0; j < kmerLength; j++) {
                tmp[j] = sequence[offset+o+j];
            }
            long bk = DiscriminativeKmers.charsToLong(tmp);
            long ak = basicKmers * ((stranded ? maxOffset : 0 ) + o - 1) + bk ;
            output[(int)ak] = 1;
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
                saveLine(training, 2, bghits.get(s), kmerFeatures(background.get(s)));
            } else {
                saveLine(test,2, bghits.get(s), kmerFeatures(background.get(s)));
                testregions.println(s);
            }
        }
        for (WeightMatrix m : matrices) {
            featurenames.println(m.getName() + "\t" + m.getVersion());
        }
        for (KmerFeature f : features) {
            featurenames.println(f.toString());
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
    public String toString() {return String.format("%s at %d",new String(kmer), offset);}
}