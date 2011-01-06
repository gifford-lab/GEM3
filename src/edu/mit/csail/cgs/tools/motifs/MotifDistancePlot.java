package edu.mit.csail.cgs.tools.motifs;

import java.io.*;
import java.util.*;
import java.sql.*;
import java.text.ParseException;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;
import edu.mit.csail.cgs.utils.database.DatabaseFactory;
import edu.mit.csail.cgs.utils.database.DatabaseException;
import edu.mit.csail.cgs.utils.database.UnknownRoleException;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.tools.utils.Args;


/**
 * Reads a list of Points or Regions on stdin.  These should be ordered and will
 * be presented along the x-axis of the plot.  The command line arguments specify a motif
 * and a maximum distance from the point (or center of region) to look for the motif.
 * The y-value of the plot is the average distance to a motif over a sliding window
 * of n regions/points.  The x-y values for the plot are produced on
 * stdout.
 *
 * java MotifDistancePlot --species "$MM;mm9" --wm "Cdx2;Jaspar" [--maxdistance 500] [--smooth 100] [--minpercent .75] < regions.txt > plot.txt
 */

public class MotifDistancePlot {

    public static void main(String args[]) throws Exception {
        int smooth = Args.parseInteger(args,"smooth",100);
        int maxdistance = Args.parseInteger(args,"maxdistance",500);
        double minpercent = Args.parseDouble(args,"minpercent",.75);
        Genome genome = Args.parseGenome(args).cdr();
        SequenceGenerator seqgen = new SequenceGenerator(genome);
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              genome.getDBID());
        Collection<WeightMatrix> matrices = Args.parseWeightMatrices(args);
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }
        
        int sumdistances = 0;
        int distances[] = new int[smooth]; // circular buffer for whether or not there was a motif under a region
        int temp[] = new int[smooth];
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        int lineno = 0;
        while ((line = reader.readLine()) != null) {
            Region r = null;
            try {
                r = Region.fromString(genome, line);
            } catch (Exception e) {
                System.err.println("Can't parse line " + line);
                continue;
            }
            if (smooth != 0 && lineno >= smooth && lineno % smooth == 0) {
                double sum = 0;
                for (int j = 0; j < temp.length; j++) {
                    temp[j] = distances[j];
                    sum += distances[j];
                }
                Arrays.sort(temp);
                
                int median = temp[temp.length / 2];
                double mean = sum / temp.length;

                System.out.println(String.format("%d\t%d\t%.2f",
                                                 (lineno - smooth),
                                                 median, mean));
            }
            if (r == null) {
                System.err.println("Null region from " + line);
                continue;
            }

            int center = (r.getStart() + r.getEnd()) / 2;
            int start = center - maxdistance;
            int end = center + maxdistance;
            if (start < 0) { start = 0; }
            if (end > genome.getChromLength(r.getChrom())) { end = genome.getChromLength(r.getChrom()) - 1;}
            r = new Region(r.getGenome(), r.getChrom(), start, end);
            char string[] = seqgen.execute(r).toCharArray();
            int mindistance = maxdistance * 2;
            for (WeightMatrix wm : matrices) {
                float minscore = (float)(wm.getMaxScore() * minpercent);
                List<WMHit> hits = WeightMatrixScanner.scanSequence(wm,
                                                                    minscore,
                                                                    string);
                for (WMHit h : hits) {
                    int c = start + (h.start + h.end) / 2;
                    int d = Math.abs(c - center);
                    if (d < mindistance) {
                        mindistance = d;
                    }
                }
            }
            if (smooth == 0) {
                System.out.println(mindistance);
                continue;
            }

            int oldval = distances[lineno % smooth];
            if (lineno >= smooth) {
                sumdistances -= oldval;
            }
            sumdistances += mindistance;
            distances[lineno % smooth] = mindistance;
            lineno++;
        }


    }


}