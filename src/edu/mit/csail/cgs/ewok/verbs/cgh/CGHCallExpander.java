package edu.mit.csail.cgs.ewok.verbs.cgh;

import java.util.*;
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

/** expands an input region (the query) into a set of CGH calls for
 *   regions that have increased or decreased in copy numbeer.  The score
 *   of each output region indicates the log of the multiplicative change from the Cy5 channel
 *   to the Cy3 channel.
 *
 *   eg, copy number 1 to 2 gives a score of 1
 *       2 to 3 gives a score of 0.4
 *
 */


public class CGHCallExpander implements Expander<Region, ScoredRegion> {

    /*states are 0: low ratio
                 1: ratio = 1
                 2: high ratio
                 
      all ratios are log-ratios
    */
    public static final int numStates = 3;
    public static final int LOW = 0, ONE = 1, HIGH = 2;
    private static final String[] stateNames = {"low","one","high"};
    /* default hmm parameters.
        mean then variance
    */
    public static final double[][] defaultParams = {{-.42, .3},
                                                    {0, .3},
                                                    {.42, .3}};

    private ChipChipData data;
    private Hmm<ObservationReal>  hmm;
    private ArrayList<Opdf<ObservationReal>> pdfs;
    private double[][] transition;
    private double[] initialProbabilities;

    public CGHCallExpander (ChipChipData data) {
        this.data = data;

        transition = new double[numStates][numStates];
        transition[LOW][LOW] = .98;
        transition[LOW][ONE] = .02;
        transition[LOW][HIGH] = .00;
        transition[ONE][LOW] = .02;
        transition[ONE][ONE] = .94;
        transition[ONE][HIGH] = .02;
        transition[HIGH][LOW] = .00;
        transition[HIGH][ONE] = .02;
        transition[HIGH][HIGH] = .98;
        initialProbabilities = new double[3];
        initialProbabilities[LOW] = .01;
        initialProbabilities[ONE] = .98;
        initialProbabilities[HIGH] = .01;

        pdfs = new ArrayList<Opdf<ObservationReal>>();
        /* params are mean and variance */
        pdfs.add(new OpdfGaussian(defaultParams[0][0],
                                  defaultParams[0][1]));
        pdfs.add(new OpdfGaussian(defaultParams[1][0],
                                  defaultParams[1][1]));
        pdfs.add(new OpdfGaussian(defaultParams[2][0],
                                  defaultParams[2][1]));
        hmm = new Hmm<ObservationReal>(initialProbabilities,
                                       transition,
                                       pdfs);                              
    }

    /* learns the mean and variance for the specified state
       from the provided regions.  This overwrites the old 
       PDF for that state
    */
    public void learnStateFromRegions(int state,
                                      Collection<Region> regions) throws NotFoundException {
        DoubleArrayList values = new DoubleArrayList();
        for (Region r : regions) {
            data.window(r.getChrom(), r.getStart(), r.getEnd());
            for (int i = 0; i < data.getCount(); i++) {
                for (int j = 0; j < data.getReplicates(i); j++) {
                    double d = Math.log(data.getRatio(i,j));
                    if (!(Double.isInfinite(d) || Double.isNaN(d))) {
                        values.add(d);
                    }
                }
            }
        }
        double mean = Descriptive.mean(values);
        double std = Descriptive.standardDeviation(Descriptive.variance(values.size(),
                                                                        Descriptive.sum(values),
                                                                        Descriptive.sumOfSquares(values)));
        /* TODO : verify that OpdfGaussian wants the variance rather than the stddev */
        pdfs.set(state, new OpdfGaussian(mean, std * std));
    }
    /* set parameters for a state.  This is the mean and standard deviation of the log ratios for that state*/
    public void setStateParameters(int state,
                                   double mean,
                                   double stddev) {
        pdfs.set(state, new OpdfGaussian(mean, stddev * stddev));
    }

    /* after calling setStateParaemters or learnStateFromRegions, you must
       call resetHMM() to actually update the internal hidden markov model
    */
    public void resetHMM() {
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
                double d = Math.log(data.getRatio(i,j));
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
            if (states[i] != ONE) {
                int j = i;
                double ratiosum = 0;
                while (j + 1 < states.length && states[j+1] == states[i]) {
                    ratiosum += observations.get(j).value;
                    j++;
                }
                output.add(new ScoredRegion(r.getGenome(),
                                            r.getChrom(),
                                            positions.get(i),
                                            positions.get(j),
                                            ratiosum / (j - i + 1)));
                i = j;
            }
        }
        return output;
    }


}