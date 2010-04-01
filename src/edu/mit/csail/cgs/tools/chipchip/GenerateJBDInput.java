package edu.mit.csail.cgs.tools.chipchip;

import java.io.*;
import java.util.*;
import java.sql.*;
import cern.colt.list.DoubleArrayList;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.utils.io.parsing.textfiles.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.tools.utils.Args;

/**
 * Generates input datafiles for JBD.
 * Usage:
 *   java edu.mit.csail.cgs.tools.chipchip.GenerateJBDInput --species "$SC;SGDv1" --expt "name;version;replicate" --outname base_filename 
 * Optional arguments
 * --minstd .3 : minimum standard deviation for a probe
 * --coeffsonly : only dump the coefficients file
 * --maxdist 700 : don't dump the fragment length distribution beyond this size
 * --unique 3 : don't dump probes that have more than this many genomic locations
 * --design 'Sc 244k' : only dump probes from this array design
 * --pieces : break the output file into pieces at large gaps in the probes.  "large gap" is 2X the fragdist length
 *
 * Output is tab separated columns of
 * - position (center of probe)
 * - cy5
 * - cy3
 * - ratio
 * - stddev across replicates
 * - replicate #
 */
public class GenerateJBDInput {

    private int maxdist, unique, neededGap;
    private double minstd;
    private boolean coeffsonly, pieces;
    private String outname, designName;
    private List<Integer> exptIDs, fragdistIDs;
    private Map<Integer,Integer> exptToColumn;
    private int designid;
    private ChipChipMetadataLoader loader;
    private Genome genome;
    private java.sql.Connection cxn;
    

    public GenerateJBDInput() throws SQLException {
        loader = new ChipChipMetadataLoader();
    }

    public void parseArgs(String args[]) throws NotFoundException, SQLException {
        genome = Args.parseGenome(args).cdr();                                 
        maxdist = Args.parseInteger(args,"maxdist",-1);
        unique = Args.parseInteger(args,"unique",-1);
        minstd = Args.parseDouble(args,"minstd",-1);
        outname = Args.parseString(args,"outname",null);
        coeffsonly = Args.parseFlags(args).contains("coeffsonly");
        pieces = Args.parseFlags(args).contains("pieces");
        designName = Args.parseString(args, "design", null);
        List<ExptNameVersion> envs = Args.parseENV(args,"expt");        
        exptToColumn = new HashMap<Integer,Integer>();
        fragdistIDs = new ArrayList<Integer>();
        exptIDs = new ArrayList<Integer>();
        cxn = DatabaseFactory.getConnection("chipchip");
        int i = 0;
        for (ExptNameVersion env : envs) {
            Collection<Experiment> expts = loader.loadExperiment(env);
            for (Experiment e : expts) {
                exptIDs.add(e.getDBID());
                fragdistIDs.add(e.getFragDist());
                exptToColumn.put(e.getDBID(), ++i);
            }
        }       
    }

    public void dumpCoeffs() throws IOException, SQLException {
        PreparedStatement ps = cxn.prepareStatement("select distance, value from fragdistentry where distribution = ? order by distance");
        neededGap = 0;
        for (int i = 0; i < fragdistIDs.size(); i++) {
            int id = fragdistIDs.get(i);
            ps.setInt(1,id);
            ResultSet rs = ps.executeQuery();
            String fname = String.format("%s.%d.coeffs", outname, i + 1);
            PrintWriter pw = new PrintWriter(new FileOutputStream(new File(fname)));
            while (rs.next()) {
                int dist = rs.getInt(1);
                double val = rs.getDouble(2);
                if (maxdist > 0 && dist > maxdist) {
                    break;
                }
                if (dist > neededGap) {
                    neededGap = dist;
                }
                pw.println(String.format("%d\t%f", dist, val));
            }
            rs.close();
            pw.close();
        }
        ps.close();
    }

    public void dumpData() throws IOException, SQLException, NotFoundException {
        String uniqueclause = unique > 0 ? (" and probelocation.loccount <= " + unique) : "";
        StringBuffer exptclause = new StringBuffer();
        exptclause.append(" data.experiment in (");
        for (int i = 0; i < exptIDs.size(); i++) {
            if (i == 0) {
                exptclause.append(exptIDs.get(i));
            } else {
                exptclause.append("," + exptIDs.get(i));
            }
        }
        exptclause.append(")");
        ResultSet rs;
        Statement stmt = cxn.createStatement();
        if (designName != null) {
            ArrayDesign design = loader.getArrayDesign(designName, genome);
            int designid = design.getDBID();
            rs = stmt.executeQuery("select probelocation.chromosome, data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, data.experiment " +
                                   " from data, probelocation, probedesign where " + exptclause.toString() + " and data.probe = probelocation.id " +
                                   " and data.probe = probedesign.id and probedesign.arraydesign =  " + designid + 
                                   uniqueclause + 
                                   " order by probelocation.chromosome, probelocation.startpos, data.experiment ");
        } else {
            String cmd = "select /*+ LEADING (data)*/ /*+ INDEX (ix_data_expt)*/ " +
                " probelocation.chromosome, data.channelone, data.channeltwo, data.ratio, probelocation.startpos, probelocation.stoppos, data.experiment " +
                " from data, probelocation where " + exptclause.toString() + " and data.probe = probelocation.id  " +
                uniqueclause +                                   
                " order by probelocation.chromosome, probelocation.startpos, data.experiment ";
            rs = stmt.executeQuery(cmd);
        }
        int lastchrom = -1;
        int lastpos = -1;
        ArrayList<String> lines = new ArrayList<String>();
        int piece = 1;
        PrintWriter file = null;

        DoubleArrayList ratios = new DoubleArrayList();
        List<Row> rows = new ArrayList<Row>();
        
        while (rs.next()) {
            Row r = new Row(rs);
            if (r.startpos != lastpos || r.chrom != lastchrom ) {
                if (file != null) {
                    double mean = cern.jet.stat.Descriptive.mean(ratios);
                    double stddev = Math.sqrt(cern.jet.stat.Descriptive.sampleVariance(ratios,mean));
                    if (stddev < minstd) {
                        stddev = minstd;
                    }
                    for (Row row : rows) {
                        file.println(String.format("%d\t%f\t%f\t%f\t%f\t%d",
                                                   (row.startpos + row.stoppos)/2,
                                                   row.cy5,
                                                   row.cy3,
                                                   row.ratio,
                                                   stddev,
                                                   exptToColumn.get(row.expt)));
                    }
                }
                ratios.clear();
                rows.clear();
            }

            if ((r.chrom != lastchrom) ||
                (pieces && (r.startpos - lastpos > 2*neededGap))) {
                System.err.println("New Chrom " + r.chrom);
                String chromname = null;
                try {
                    chromname = genome.getChromName(r.chrom);
                } catch (NullPointerException e) {
                    /* this means the chromosome is from another genome.  That's fine;
                       array data gets mapped to multiple genomes.  Just skip it
                    */
                }

                if (chromname != null) {
                    if (r.chrom != lastchrom) {
                        piece = 1;
                    } else {
                        piece++;
                    }
                    String fname = String.format("%s.%s.%s.data",
                                                 outname, chromname,piece);
                    if (file != null) {
                        file.close();
                    }
                    System.err.println("NEW FILE " + fname);
                    file = new PrintWriter(new FileOutputStream(new File(fname)));
                } else {
                    file = null;
                }

            }

            rows.add(r);
            ratios.add(r.ratio);
            lastchrom = r.chrom;
            lastpos = r.startpos;

        }
        if (file != null) {
            file.close();        
        }
    }

    public void dumpArgs() throws IOException, NotFoundException, SQLException {
        PrintWriter pw = new PrintWriter(new FileOutputStream(new File("./args.txt")));
        pw.print(String.format("--species '%s;%s' --pattern '%s' ",
                               genome.getSpecies(),
                               genome.getVersion(),
                               outname));
        for (int id : exptIDs) {
            Experiment e = loader.loadExperiment(id);
            pw.print(String.format(" --expt '%s;%s;%s' ",
                                   e.getName(), e.getVersion(), e.getReplicate()));
        }
        
        pw.close();
    }

    public static void main(String[] args) throws Exception {
        GenerateJBDInput g = new GenerateJBDInput();
        g.parseArgs(args);
        g.dumpCoeffs();
        if (!g.coeffsonly) {
            g.dumpData();
        }
        g.dumpArgs();
    }
}

class Row {
    public int chrom, startpos, stoppos, expt;
    public double cy5, cy3, ratio;
    public Row (ResultSet rs) throws SQLException {
        chrom = rs.getInt(1);
        cy5 = rs.getDouble(2);
        cy3 = rs.getDouble(3);
        ratio = rs.getDouble(4);
        startpos = rs.getInt(5);
        stoppos = rs.getInt(6);
        expt = rs.getInt(7);
    }
}