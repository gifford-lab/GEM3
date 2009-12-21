package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import java.sql.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.alignments.MultiZAlignRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.utils.database.*;

/**
 * Maps a region to a set of regions that describe the alignment of
 * the input region to the specified genome.  The
 * output regions are regions in the other genome.  These alignments
 * are in the UCSC annotation tables (or use the same format).
 *
 * MultiZAlign generator can use different sets of alignments as
 * input.  The default set is chain, where the tablename is
 * <chromosome>_chain<othergenome>.  getAlignmentVerions will return a
 * set of all alignments that can be used for a set of genomes, based
 * on the edu.mit.csail.cgs.ewok.species_alignments file.
 *
 * The aligned regions can be a little confusing.  The input to
 * execute() is matched against the TARGET region in the alignments.
 * In the resulting MultiZAlignRegion, the base
 * genome/chrom/start/stop (the values you get back with getGenome,
 * getChrom, etc) are the TARGETs from the alignment and the other
 * genome/chrom/start/stop (the values you get back with
 * getOtherGenome, getOtherChrom) are the QUERY from the alignment.
 */

public class MultiZAlignGenerator<X extends Region> implements Expander<X,MultiZAlignRegion> {

    public static String defaultAlignmentPrefix = "chain";
    private Genome genome, other;
    private String otherGenome, alignPrefix;

    public MultiZAlignGenerator(Genome genome, Genome other) {
        this.genome = genome;
        this.other = other;
        this.otherGenome = 
            Character.toUpperCase(other.getVersion().charAt(0)) + 
            other.getVersion().substring(1);
        this.otherGenome = this.otherGenome.replaceAll("[^\\w]","_");
        alignPrefix = defaultAlignmentPrefix;
    }

    /* The alignPrefix lets a MultiZAlignGenerator retrieve data from different sets of tables
       during runtime.  This allows different data sources, ie different sets of alignments between
       genomes */
    public String getAlignPrefix() {return alignPrefix;}
    public void setAlignPrefix(String p) {alignPrefix = p;}

    public Iterator<MultiZAlignRegion> execute(X region) {
        String chr = region.getChrom();
        if (!chr.matches("^(chr|scaffold).*")) {
            chr = "chr" + chr;
        }
        String tablename = chr + "_" + alignPrefix + otherGenome;
        try {
            java.sql.Connection cxn =
                genome.getUcscConnection();
            PreparedStatement ps = cxn.prepareStatement("select replace(qName,'chr',''), qStart, qEnd, replace(tName,'chr',''), tStart, tEnd, score, qStrand from " + tablename + " where " +
                                                        "((tStart <= ? and tEnd >= ?) or (tStart >= ? and tStart <= ?)) order by tStart");            
            ps.setInt(1,region.getStart());
            ps.setInt(2,region.getStart());
            ps.setInt(3,region.getStart());
            ps.setInt(4,region.getEnd());
            ResultSet rs = ps.executeQuery();
            ArrayList<MultiZAlignRegion> results = new ArrayList<MultiZAlignRegion>();
            while (rs.next()) {
                MultiZAlignRegion r = new MultiZAlignRegion(genome,
                                                            rs.getString(4),
                                                            rs.getInt(5),
                                                            rs.getInt(6),
                                                            other,
                                                            rs.getString(1),
                                                            rs.getInt(2),
                                                            rs.getInt(3),
                                                            rs.getDouble(7),
                                                            rs.getString(8).charAt(0));
                results.add(r);
            }
            rs.close();
            ps.close();
            DatabaseFactory.freeConnection(cxn);
            return results.iterator();
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't get UCSC MultiZ for " + tablename,ex);
        }
    }

    static Map<String,List<String>> alignmentVersions = null;
    public static Map<String,List<String>> getAlignmentVersions() {
        if (alignmentVersions != null) {
            return alignmentVersions;
        }
        ResourceBundle res = ResourceBundle.getBundle("edu.mit.csail.cgs.ewok.species_alignments");
        Map<String,List<String>> map = new HashMap<String,List<String>>();
        Enumeration<String> keys = res.getKeys();
        while (keys.hasMoreElements()) {
            String k = keys.nextElement();
            String[] pieces = k.split("___");
            if (pieces.length < 3) {
                System.err.println("TOO SHORT : " + k);
                continue;
            }
            String genomes = pieces[0] + "\t" + pieces[1];
            List<String> methods;
            if (!map.containsKey(genomes)) {
                methods = new ArrayList<String>();
                map.put(genomes,methods);
            } else {
                methods = map.get(genomes);
            }
            methods.add(pieces[2]);            
        }
        alignmentVersions = map;
        return map;
    }
    public static Set<String> getAlignmentVersions(Set<String> genomes) {
        Map<String,List<String>> all = getAlignmentVersions();
        for (String k : all.keySet()) {
            System.err.println("Found pair:" + k);
        }

        Map<String,Integer> counts = new HashMap<String,Integer>();
        for (String g1 : genomes) {
            for (String g2 : genomes) {
                String k = g1.replaceAll("[^\\w_]+","_") + "\t" + g2.replaceAll("[^\\w_]+","_");
                if (g1.equals(g2)) {continue;}
                if (!all.containsKey(k)) {
                    System.err.println("No Entry for " + k);
                    continue;
                }
                for (String a : all.get(k)) {
                    if (!counts.containsKey(a)) {
                        counts.put(a,1);
                    } else {
                        counts.put(a,counts.get(a) + 1);
                    }
                }
            }
        }
        Set<String> output = new HashSet<String>();
        int target = genomes.size() * (genomes.size() - 1);
        for (String a : counts.keySet()) {
            if (counts.get(a) == target) {
                output.add(a);
            }
        }
        return output;
    }
}
