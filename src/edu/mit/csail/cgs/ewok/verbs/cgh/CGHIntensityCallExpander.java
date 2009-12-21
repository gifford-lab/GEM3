package edu.mit.csail.cgs.ewok.verbs.cgh;

import java.util.*;
import java.sql.*;
import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.database.*;
import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;
import edu.mit.csail.cgs.datasets.chipchip.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * Performs an intensity-based single channel analysis of CGH data.  Starts
 * by computing the median intensity and uses this as a baseline intensity.
 * Generates an HMM for 0-copies (intensity lower than baseline), 1-copy
 * (intensity is baseline), 2 copies (2x baseline), etc
 *
 * The output are the regions in which the count != 1.  The score for the output
 * regions is average copy number over the region (this is going to be the average
 * intensity in the ouptut region divided by the median intensity in the in the entire
 * experiment).
 */

public class CGHIntensityCallExpander implements Expander<Region,ScoredRegion> {

    private ChipChipData data;
    private Hmm<ObservationReal>  hmm;
    private ArrayList<Opdf<ObservationReal>> pdfs;
    private double[][] transition;
    private double[] initialProbabilities;
    private double medianIntensity;
    private int numStates; // number of states in the HMM
    private double stdConst; // intensity * stdconst = stddev of the observation
    private boolean cy5;
    private double maxIntensity; // highest intensity in the HMM; don't let observations go above this or you may get numeric overflow in the HMM
    private double selfTransition = .99999999999;
    private boolean computeOnLog = true;

    /** Create a new CGHCallExpander.  The experiment id provided
     * must correspond to the ChipChipData provided
     */
    public CGHIntensityCallExpander (int exptid, ChipChipData data, boolean cy5) {
        this.data = data;
        this.cy5 = cy5;
        init(exptid, data);
    }
    /** Create a new CGHCallExpander.  The experiment id provided
     * must correspond to the ChipChipData provided.  selfTransitionProbability
     * is the probability of transitioning from a state to itself in the HMM.
     * The default value is .99999999999999; values closer to one should make the
     * HMM more resistant to single-probe noise.
     */
    public CGHIntensityCallExpander (int exptid, ChipChipData data, boolean cy5, double selfTransitionProbability) {
        this.data = data;
        this.cy5 = cy5;
        this.selfTransition = selfTransitionProbability;
        init(exptid, data);
    }

    private void init(int exptid, ChipChipData data) {
        try {
            java.sql.Connection chipcxn = DatabaseFactory.getConnection("chipchip");
            PreparedStatement stmt = chipcxn.prepareStatement("select median(channelone), median(channeltwo) from data where experiment = ?");
            stmt.setInt(1,exptid);
            ResultSet rs = stmt.executeQuery();
            if (rs.next()) {
                medianIntensity = rs.getDouble(cy5 ? 1 : 2);
            } else {
                throw new DatabaseException("no next");                    
            }
        } catch (UnknownRoleException ex) {
            throw new DatabaseException("Couldn't connect with role chipchip",ex);
        } catch (SQLException ex) {
            throw new DatabaseException("Couldn't connect to database for role chipchip",ex);
        }
        System.err.println("Learned median intensity as " + medianIntensity);
        
        stdConst = .3;
        numStates = 6;
        initialProbabilities = new double[numStates];
        transition = new double[numStates][numStates];
        pdfs = new ArrayList<Opdf<ObservationReal>>();
        double means[] = new double[numStates];
        for (int i = 0; i < numStates; i++) {
            means[i] = computeOnLog ? 
                Math.log(Math.max(medianIntensity * i, 1000)) :
                Math.max(medianIntensity * i, 1000);
        }
        for (int i = 0; i < numStates; i++) {
            initialProbabilities[i] = (i == 1) ? .9 : (.1 / (numStates - 1));
            double std;
            if (computeOnLog) {
                if (i == 0) {
                    std = (means[1] - means[0]) / 4;
                } else if (i == numStates - 1) {
                    std = (means[numStates - 1] - means[numStates - 2]) / 4;
                } else {
                    std = Math.min(means[i] - means[i-1],
                                   means[i+1] - means[i]) / 4;
                }
            } else {
                std = means[i] * stdConst;
            }
            System.err.println(String.format("state %d, mean %.2f, stddev %.2f",
                                             i, means[i], std));
            pdfs.add(new OpdfGaussian(means[i],std * std));

            for (int j = 0; j < numStates; j++) {
                if (i == j) {
                    transition[i][j] = selfTransition;
                } else {
                    transition[i][j] = (1 - selfTransition) / (numStates - 1);
                }
            }
        }
        maxIntensity = medianIntensity * (numStates + 2);
        hmm = new Hmm<ObservationReal>(initialProbabilities,
                                       transition,
                                       pdfs);                              
    }

    public Pair<List<ObservationReal>, List<Integer>> getObservations(Region r) throws NotFoundException {
        data.window(r.getChrom(),
                    r.getStart(),
                    r.getEnd());
        ArrayList<ObservationReal> observations = new ArrayList<ObservationReal>();
        ArrayList<Integer> positions = new ArrayList<Integer>();
        for (int i = 0; i < data.getCount(); i++) {
            for (int j = 0; j < data.getReplicates(i); j++) {
                double d = Math.min(cy5 ? data.getIP(i,j) : data.getWCE(i,j), maxIntensity);
                if (computeOnLog) {
                    d = Math.log(d);
                }
                if (!(Double.isInfinite(d) || Double.isNaN(d))) {
                    observations.add(new ObservationReal(d));
                    positions.add(data.getPos(i));
                }
            }
        }
        return new Pair<List<ObservationReal>, List<Integer>>(observations, positions);
    }

    public Iterator<ScoredRegion> execute(Region r)  {
        return executeList(r).iterator();
    }
    public List<ScoredRegion> executeList(Region r) {
        Pair<List<ObservationReal>, List<Integer>> pair;
        try {
            pair = getObservations(r);
        } catch (NotFoundException e) {
            throw new DatabaseException(e.toString(), e);
        }
        List<ObservationReal> observations = pair.car();
        List<Integer> positions = pair.cdr();
        ViterbiCalculator vc = new ViterbiCalculator(observations,hmm);
        int[] states = vc.stateSequence();
        List<ScoredRegion> output = new ArrayList<ScoredRegion>();
        for (int i = 0; i < states.length; i++) {
//             System.err.println(String.format("%d : %d --> %d, %s",
//                                              i, positions.get(i),
//                                              states[i],
//                                              observations.get(i).toString()));
            if (states[i] != 1) {
                int j = i;
                double intensitysum = observations.get(i).value;               
                while (j + 1 < states.length && states[j+1] == states[i]) {
//                     System.err.println(String.format("%d : %d --> %d, %s",
//                                                      j, positions.get(j),
//                                                      states[j],
//                                                      observations.get(j).toString()));
                    intensitysum += observations.get(j).value;
                    j++;
                }
                output.add(new ScoredRegion(r.getGenome(),
                                            r.getChrom(),
                                            positions.get(i),
                                            positions.get(j),
                                            (intensitysum / (j - i + 1)) / medianIntensity));
                i = j;
            }
        }
        return output;
    }       

}