package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Color;
import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import edu.mit.csail.cgs.datasets.motifs.*;

import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.ewok.verbs.*;


public class DiscriminativeKmers {

    private Genome genome;
    private double minfoldchange;
    private int k, mask, maxmismatch, minclustersize, minclustercount;
    private Binomial binomial;
    private SequenceGenerator seqgen;
    private int randombgcount = 1000,  // number of random background regions to pic
        randombgsize = 100, // size of random background regions
        parsedregionexpand; // expand input regions by this much on either side
    private Map<String,char[]> foreground, background;  // foreground and background sequences
    private boolean printKmers;
    private List<WeightMatrix> pwms;
    private String outbase;
    private final static long intmask = 0xffffffffL;
    private final static int maxshift = 3;
    private final static char[] toChar = {'A','C','G','T'};

    public static long charsToLong(char[] chars) {
        long out = 0;
        for (int i = 0; i < chars.length; i++) {
            out <<= 2;
            char newchar = chars[i];
            switch (newchar) {
            case 'A':
            case 'a':
                out += 0;
            break;
            case 'C':
            case 'c':
                out += 1;
            break;
            case 'G':
            case 'g':
                out += 2;
            break;
            case 'T':
            case 't':
                out += 3;
            break;

            default:
                break;
            }
        }
        return out;
    }
    /** adds another character to the end (right side) of the long representation
     * of a kmer.  mask is the bitmask of usable bits in the long rep.
     */
    public static long addChar(long existing, long mask, char newchar) {
        existing <<= 2;
        switch (newchar) {
        case 'A':
        case 'a':
            existing += 0;
        break;
        case 'C':
        case 'c':
            existing += 1;
        break;
        case 'G':
        case 'g':
            existing += 2;
        break;
        case 'T':
        case 't':
            existing += 3;
        break;

        default:
            break;
        }
        return existing & mask;
    }
    /** converts a long representation of a kmer back to a string. 
     */
    public static String longToString(long l, int k) {
        return new String(longToChars(l,k));
    }
    public static  char[] longToChars(long l, int k) {
        char[] chars = new char[k];
        while (k-- > 0) {
            int b = (int)(l & 3);
            l >>= 2;
            chars[k] = toChar[b];
        }
        return chars;
    }
    public static long reverseComplement(long kmer,
                                         int k) {
        long out = 0;
        for (int i = 0; i < k; i++) {
            byte b = (byte)((kmer ^ 3) & 3);
            kmer >>= 2;
            out = (out << 2) | b;
        }
        return out;
    }
    /** returns a map from kmers to their frequency 
     * in the string.
     */
    public static Map<Long,Integer> count(char[] chars, 
                                          int k,
                                          long mask,
                                          Map<Long,Integer> map) {
        if (map == null) {
            map = new HashMap<Long,Integer>();
        }
        long l = 0;
        for (int i = 0; i < k-1; i++) {
            l = addChar(l,mask,chars[i]);
        }
        for (int i = k; i < chars.length; i++) {
            l = addChar(l,mask,chars[i]);
            if (!map.containsKey(l)) {
                map.put(l,1);
            } else {
                map.put(l,map.get(l)+1);
            }            
        }
        return map;
    }
    /** return the number of bases that two kmers have in common.
     */
    public static short countBasesSame(long a, long b, int k) {
        long x = a ^ b;
        short same = 0;
        while (k-- > 0) {
            if ((x & 3) == 0) {
                same++;
            }
            x >>= 2;
        }
        return same;
    }
    /** returns a long describing the best match between a and b.  The
     * upper int is the shift- the offset at which a best matched b.
     * Positive shift means to shift a to the right such that the last
     * characters of a and first characters of b are overhanging.The
     * lower int is the number of bases in common.
     */
    public static int countBasesSameOneDir(long a, long b, int k, int maxShift) {
        short bestSame = countBasesSame(a,b,k);
        short bestPos = 0;
        for (short l = 1; l <= maxShift; l++) {
            long newa = a >> l*2;
            short s = countBasesSame(newa,b,k-l);
            if (s > bestSame) { 
                bestSame = s;
                bestPos = l;
            }
        }
        for (short r = 1; r <= maxShift; r++) {
            long newb = b >> r*2;
            short s = countBasesSame(a,newb,k-r);
            if (s > bestSame) { 
                bestSame = s;
                bestPos = (short)(-1 * r);
            }
        }
        int ret = bestPos;
        ret = (ret << 16) | bestSame;
        return ret;
    }
    /** first short is zero
        second short is one if reverse-complement was best
        third short is offset
        fourth short (LSB) is score
    */
    public static long countBasesSame(long a, long b, int k, int maxShift) {
        long forw = countBasesSameOneDir(a,b,k,maxShift);
        long rc = countBasesSameOneDir(reverseComplement(a,k),b,k,maxShift);
        if ((forw & 0xffff) > (rc & 0xffff)) {
            return (forw << 1);
        } else {
            return (rc << 1) | 1;
        }
    }
    public static short getRC(long l) {
        return (short)(l & 1);
    }
    public static short getSameness(long l) {
        return (short)((l >> 1) & 0xffff);
    }
    public static short getShift(long l) {
        return (short)((l >> 17) & 0xffff);
    }
    public static void paintMotif(WeightMatrix wm, String fname) throws IOException {
        File f = new File(fname);
        BufferedImage im = 
            new BufferedImage(800, 200, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2 = im.createGraphics();
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        WeightMatrixPainter wmp = new WeightMatrixPainter();
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,800,200);
        
        wmp.paint(wm,g2,0,0,800,200);
        ImageIO.write(im, "png", f);
    }
    public WeightMatrix toWeightMatrix(int sums[], int counts[][]) {
        double[] bgmodel = new double[4];
        for (String s : background.keySet()) {
            char[] chars = background.get(s);
            for (int i = 0; i < chars.length; i++) {
                switch (chars[i]) {
                case 'A':
                case 'a':
                    bgmodel[0]++;
                break;
                case 'C':
                case 'c':
                    bgmodel[1]++;
                break;
                case 'G':
                case 'g':
                    bgmodel[2]++;
                break;
                case 'T':
                case 't':
                    bgmodel[3]++;
                break;
                
                default:
                    break;
                }
            }
        }
        double bgsum = bgmodel[0] + bgmodel[1] + bgmodel[2] + bgmodel[3];
        for (int i = 0; i < bgmodel.length; i++) {
            bgmodel[i] = bgmodel[i] / bgsum;
        }
        WeightMatrix out = new WeightMatrix(sums.length);
        for (int i = 0; i < sums.length; i++) {
            out.matrix[i]['A'] = (float)Math.log( ((double)counts[i][0] / (double)sums[i]) / bgmodel[0] );
            out.matrix[i]['C'] = (float)Math.log( ((double)counts[i][1] / (double)sums[i]) / bgmodel[1] );
            out.matrix[i]['G'] = (float)Math.log( ((double)counts[i][2] / (double)sums[i]) / bgmodel[2] );
            out.matrix[i]['T'] = (float)Math.log( ((double)counts[i][3] / (double)sums[i]) / bgmodel[3] );
        }
        return out;
    }
    
    public List<KmerCluster> cluster(List<KmerCount> kmers, int k, int maxmismatch) {
        Collections.sort(kmers, new KmerCountComparator());
        List<KmerCluster> clusters = new ArrayList<KmerCluster>();
        clusters.add(new KmerCluster(kmers.remove(0)));
        for (int i = 0; i < kmers.size(); i++) {
            KmerCount kc = kmers.get(i);
            int bestCluster = -1, bestSameness = -1;
            for (int j = 0; j < clusters.size(); j++) {
                KmerCluster cluster = clusters.get(j);
                long cbs = countBasesSame(kc.kmer, cluster.centroid, k, maxshift);
                short sameness = getSameness(cbs);
                if (sameness > bestSameness) {
                    bestSameness = sameness;
                    bestCluster = j;
                }
            }
            if (bestSameness > k - maxmismatch) {
                clusters.get(bestCluster).members.add(kc);
            } else {
                clusters.add(new KmerCluster(kc));
            }
        }
        Collections.sort(clusters, new KmerClusterComparator());
        return clusters;
    }
    public void printClusters (List<KmerCluster> clusters) {
        int clusterNumber = 0;
        pwms = new ArrayList<WeightMatrix>();
        for (KmerCluster cluster : clusters) {
            if (cluster.totalCount() < minclustercount) {
                continue;
            }
            if (cluster.members.size() < minclustersize) {
                continue;
            }

            int sums[] = new int[k + 2*maxshift];
            int counts[][] = new int[k + 2*maxshift][4];
            for (int j = 0; j < sums.length; j++) {
                sums[j] = 4;
                for (int i = 0; i < 4; i++) {
                    counts[j][i] = 1;
                }
            }
            for (KmerCount kc : cluster.members) {
                long cbs = countBasesSame(kc.kmer, cluster.centroid, k, maxshift);
                short sameness = getSameness(cbs);
                short shift = getShift(cbs);
                short rc = getRC(cbs);
                long kmer = kc.kmer;
                if (rc == 1) {
                    kmer = reverseComplement(kmer,k);
                    //                    shift = (short)(shift * -1);
                }
                String kmerString = longToString(kmer,k);
                int firstbase  = maxshift + shift;
                for (int i = 1; i <= k; i++) {
                    int pos = firstbase + k - i;
                    sums[pos] += kc.count;
                    counts[pos][(int)(kmer & 3L)] += kc.count;
                    kmer >>= 2;
                }
                if (printKmers) {
                    StringBuilder sb = new StringBuilder("  ");
                    for (int i = 0; i < firstbase; i++) {
                        sb.append(" ");
                    }
                    sb.append(kmerString);
                    for (int i = 0; i < (2*maxshift-firstbase); i++) {
                        sb.append(" ");
                    }
                    sb.append("\t" + kc.count);
                    System.out.println(sb.toString());
                }
            }
            WeightMatrix wm = toWeightMatrix(sums, counts);
            wm.name = outbase + "_" + clusterNumber;
            wm.version = String.format("mfc %.2f expand %d size %d count %d",
                                       minfoldchange, parsedregionexpand, minclustersize, minclustercount);
            wm.type = "DiscriminativeKmers";

            CEResult scanresult = CompareEnrichment.doScan(wm,
                                                          foreground,
                                                          background,
                                                          null,null,
                                                          .1, .01, 2, .05, 1, null, null, false);
            if (scanresult.freqone > .1) {
                System.out.println("Cluster centroid " + longToString(cluster.centroid,k) + " and total count " + cluster.totalCount());
                for (int i = 0; i < 4; i++) {
                    System.out.print(toChar[i]);
                    for (int j = 0; j < sums.length; j++) {
                        double s = sums[j];
                        double c = counts[j][i];
                        System.out.print(String.format("\t%.2f",c/s));
                    }
                    System.out.println();
                }
                System.out.println(scanresult.toString());
                try {
                    paintMotif(wm, outbase + clusterNumber++ + ".png");
                } catch (IOException e) {
                    e.printStackTrace();
                }               
            }
        }
    }
    public DiscriminativeKmers () {
        seqgen = new SequenceGenerator();
        binomial = new Binomial(100, .01, RandomEngine.makeDefault());
    }
    public void setK(int k) {
        this.k = k;
        mask = 0;
        for (int i = 0; i < k; i++) {
            mask <<= 2;
            mask = mask | 3;
        }
    }
    public void parseArgs(String args[]) throws NotFoundException, IOException, FileNotFoundException {
        int k = Args.parseInteger(args,"k",10);
        setK(k);
        printKmers = Args.parseFlags(args).contains("printkmers");
        minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        parsedregionexpand = Args.parseInteger(args,"expand",30);
        randombgcount = Args.parseInteger(args,"randombgcount",1000);
        randombgsize = Args.parseInteger(args,"randombgsize",100);
        maxmismatch = Args.parseInteger(args,"maxmismatch",3);
        minclustersize = Args.parseInteger(args,"minclustersize",2);
        minclustercount = Args.parseInteger(args,"minclustercount",30);
        outbase = Args.parseString(args,"outbase","motif");
        genome = Args.parseGenome(args).cdr();
        String firstfname = null, secondfname = null;
        firstfname = Args.parseString(args,"first",null);
        secondfname = Args.parseString(args,"second",null);
        if (firstfname == null) {
            System.err.println("No --first specified.  Reading from stdin");
            foreground = CompareEnrichment.readRegions(genome,
                                                       new BufferedReader(new InputStreamReader(System.in)), 
                                                       parsedregionexpand,
                                                       null);
        } else {
            if (firstfname.matches(".*\\.fasta") ||
                firstfname.matches(".*\\.fa")) {
                foreground = CompareEnrichment.readFasta(new BufferedReader(new FileReader(firstfname)));
            } else {
                foreground = CompareEnrichment.readRegions(genome,
                                                           new BufferedReader(new FileReader(firstfname)), 
                                                           parsedregionexpand,
                                                           null);
            }
        }
        if (secondfname == null) {
            System.err.println("No background file given.  Generating " + randombgcount + " regions of size " + randombgsize);
            background = CompareEnrichment.randomRegions(genome,
                                                         randombgcount,
                                                         randombgsize);
        } else {
            if (secondfname.matches(".*\\.fasta") ||
                secondfname.matches(".*\\.fa")){
                background = CompareEnrichment.readFasta(new BufferedReader(new FileReader(secondfname)));
            } else {
                background = CompareEnrichment.readRegions(genome, new BufferedReader(new FileReader(secondfname)), parsedregionexpand,null);
            }
        }
    }
    public void run() {
        Map<Long,Integer> fgcounts = new HashMap<Long,Integer>();
        Map<Long,Integer> bgcounts = new HashMap<Long,Integer>();
        for (char[] chars : foreground.values()) {
            count(chars, k, mask, fgcounts);
        }
        for (char[] chars : background.values()) {
            count(chars, k, mask, bgcounts);
        }
        int fgsum = 0;
        int bgsum = 0;
        for (long l : fgcounts.keySet()) {
            fgsum += fgcounts.get(l);
        }
        for (long l : bgcounts.keySet()) {
            bgsum += bgcounts.get(l);
        }
        System.err.println("Read " + fgsum + " kmers from the fg set and " + bgsum + " from the background set");
        List<KmerCount> enriched = new ArrayList<KmerCount>();
        for (long l : fgcounts.keySet()) {
            int bgcount = (bgcounts.containsKey(l) ? bgcounts.get(l) : 0) + 1;
            double bgprob = ((double)bgcount) / ((double)bgsum);
            int fgcount = fgcounts.get(l);
            double fgprob = ((double)fgcount) / ((double)fgsum);
            
            binomial.setNandP(fgsum, bgprob);
            double pval = Math.log(1 - binomial.cdf(fgcount));

            if (fgprob > bgprob * minfoldchange) {
                KmerCount kc = new KmerCount(l,fgcount);
                // before clustering, change count from absolute counts to count above background
                kc.count = (int)(kc.count - bgprob * fgsum);
                if (kc.count > 0) {
                    enriched.add(kc);
                }
            }
        }
        List<KmerCluster> clusters = cluster(enriched, k, maxmismatch);
        printClusters(clusters);
    }
    public static void main(String args[]) throws Exception {
        DiscriminativeKmers kmers = new DiscriminativeKmers();
        kmers.parseArgs(args);
        kmers.run();
    }
}
class KmerCount {
    public long kmer;
    public int count;
    public KmerCount(long l, int c) {kmer = l; count = c;}
}
class KmerCluster {
    public long centroid;
    public List<KmerCount> members;
    public KmerCluster(KmerCount kc) {
        centroid = kc.kmer;
        members = new ArrayList<KmerCount>();
        members.add(kc);
    }
    public int totalCount() {
        int c = 0;
        for (KmerCount kc : members) {
            c += kc.count;
        }
        return c;
    }
}
class KmerCountComparator implements Comparator<KmerCount> {
    public int compare(KmerCount a, KmerCount b) {
        return b.count - a.count;
    }
}
class KmerClusterComparator implements Comparator<KmerCluster> {
    public int compare(KmerCluster a, KmerCluster b) {
        return b.totalCount() - a.totalCount();
    }
}