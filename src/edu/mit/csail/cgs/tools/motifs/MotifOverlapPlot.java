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
 * and a distance from the point (or center of region) within which an occurence of the
 * motif must fall.  The y-value of the plot is the fraction of regions near a motif; the fraction is computed
 * over a sliding window of n regions/points.  The x-y values for the plot are produced on
 * stdout.
 *
 * java MotifOverlapPlot --species "$MM;mm9" --wm "Cdx2;Jaspar" --distance 50 --smooth 100 < regions.txt > plot.txt
 */

public class MotifOverlapPlot {

    public static void main(String args[]) throws Exception {
        int smooth = Args.parseInteger(args,"smooth",100);
        int distance = Args.parseInteger(args,"distance",50);
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

        int numFound = 0;
        int found[] = new int[smooth]; // circular buffer for whether or not there was a motif under a region
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
                System.out.println(String.format("%d\t%.2f",
                                                 (lineno - smooth),
                                                 (numFound / (double)smooth)));
            }
            if (r == null) {
                System.err.println("Null region from " + line);
                continue;
            }

            int center = (r.getStart() + r.getEnd()) / 2;
            int start = center - distance;
            int end = center + distance;
            if (start < 0) { start = 0; }
            if (end > genome.getChromLength(r.getChrom())) { end = genome.getChromLength(r.getChrom()) - 1;}
            r = new Region(r.getGenome(), r.getChrom(), start, end);
            char string[] = seqgen.execute(r).toCharArray();
            boolean anyhits = false;
            for (WeightMatrix wm : matrices) {
                float minscore = (float)(wm.getMaxScore() * minpercent);
                WMHit hit = WeightMatrixScanner.scanSequenceBestHit(wm,
                                                                    minscore,
                                                                    string);
                if (hit != null) {
                    anyhits = true;
                    break;
                }
            }
            int newval = anyhits ? 1 : 0;
            if (smooth == 0) {
                System.out.println(newval);
                continue;
            }

            int oldval = found[lineno % smooth];
 
            if (oldval != 0 && lineno >= smooth) {
                numFound--;
            }
            if (newval != 0) {
                numFound++;
            }
            found[lineno % smooth] = newval;
            lineno++;
        }


    }


}