package edu.mit.csail.cgs.tools.motifs;

import java.util.*;
import java.io.*;
import java.sql.*;
import cern.jet.random.Binomial;
import cern.jet.random.engine.RandomEngine;

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
    private int k, mask;
    private Binomial binomial;
    private SequenceGenerator seqgen;
    private int randombgcount = 1000,  // number of random background regions to pic
        randombgsize = 100, // size of random background regions
        parsedregionexpand; // expand input regions by this much on either side
    private Map<String,char[]> foreground, background;  // foreground and background sequences
    private boolean printKmers, cluster;

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
    public static String longToString(long l, int len) {
        StringBuilder builder = new StringBuilder();
        while (len-- > 0) {
            int b = (int)(l & 3);
            l >>= 2;
            switch (b) {
            case 0:
                builder.append("A");
                break;
            case 1:
                builder.append("C");
                break;
            case 2:
                builder.append("G");
                break;
            case 3:
                builder.append("T");
                break;
            default:
                break;
            }
        }
        builder.reverse();
        return builder.toString();
    }
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
    public static int countBasesSame(long a, long b, int k) {
        long x = a ^ b;
        int same = 0;
        while (k-- > 0) {
            if ((x & 3) == 0) {
                same++;
            }
            x >>= 2;
        }
        return same;
    }
    public static int countBasesSame(long a, long b, int k, int maxShift) {
        int bestSame = countBasesSame(a,b,k);
        for (int l = 1; l <= maxShift; l++) {
            long newa = a >> l*2;
            int s = countBasesSame(newa,b,k-l);
            if (s > bestSame) { bestSame = s;}
        }
        for (int r = 1; r <= maxShift; r++) {
            long newb = b >> r*2;
            int s = countBasesSame(a,newb,k-r);
            if (s > bestSame) { bestSame = s;}
        }
        return bestSame;
    }
    public static void cluster(List<KmerCount> kmers, int k) {
        Collections.sort(kmers, new KmerCountComparator());
        List<KmerCluster> clusters = new ArrayList<KmerCluster>();
        clusters.add(new KmerCluster(kmers.remove(0)));
        for (int i = 0; i < kmers.size(); i++) {
            KmerCount kc = kmers.get(i);
            int bestCluster = -1, bestSameness = -1;
            for (int j = 0; j < clusters.size(); j++) {
                KmerCluster cluster = clusters.get(j);
                int sameness = countBasesSame(kc.kmer, cluster.centroid, k, 3);
                if (sameness > bestSameness) {
                    bestSameness = sameness;
                    bestCluster = j;
                }
            }
            if (bestSameness > k - 3) {
                clusters.get(bestCluster).members.add(kc);
            } else {
                clusters.add(new KmerCluster(kc));
            }
        }
        for (KmerCluster cluster : clusters) {
            System.out.println("Cluster centroid " + longToString(cluster.centroid,k) + " and total count " + cluster.totalCount());
            for (KmerCount kc : cluster.members) {
                System.out.println("  " + longToString(kc.kmer,k) + "\t" + kc.count);
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
            mask <<= 1;
            mask = mask | 1;
        }
    }
    public void parseArgs(String args[]) throws NotFoundException, IOException, FileNotFoundException {
        int k = Args.parseInteger(args,"k",10);
        setK(k);
        printKmers = Args.parseFlags(args).contains("printkmers");
        cluster = Args.parseFlags(args).contains("cluster");
        minfoldchange = Args.parseDouble(args,"minfoldchange",1);
        parsedregionexpand = Args.parseInteger(args,"expand",30);
        randombgcount = Args.parseInteger(args,"randombgcount",1000);
        randombgsize = Args.parseInteger(args,"randombgsize",100);
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
                enriched.add(new KmerCount(l,fgcount));
                if (printKmers) {
                    System.out.println(String.format("%s\t%.2e\t%.2e\t%d\t%d\t%2e\t%d\t%d\t%.2e",
                                                     longToString(l,k),
                                                     fgprob / bgprob,
                                                     pval,
                                                     fgcount,
                                                     fgsum,
                                                     fgprob,
                                                     bgcount,
                                                     bgsum,
                                                     bgprob));
                }
            }
        }
        if (cluster) {
            cluster(enriched, k);
        }
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
        return a.count - b.count;
    }
}
