package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class WeightMatrixExpander implements Expander<Region,WeightMatrixHit> {

    private Vector<WeightMatrix> matrices;
    private Vector<Float> cutoffs;
    private SequenceGenerator seqgen;
    
    public WeightMatrixExpander(WeightMatrix m, float c) {
        matrices = new Vector<WeightMatrix>();
        cutoffs = new Vector<Float>();
        if (m == null) {
            throw new NullPointerException("don't give me a null matrix!");
        }

        matrices.add(m);
        cutoffs.add(c);
        seqgen = new SequenceGenerator();
    }
    
    public WeightMatrixExpander() {
        matrices = new Vector<WeightMatrix>();
        cutoffs = new Vector<Float>();
        seqgen = new SequenceGenerator();
    }
    
    public void addWeightMatrix(WeightMatrix m, float c) { 
        if (m == null) {
            throw new NullPointerException("don't give me a null matrix!");
        }
        matrices.add(m);
        cutoffs.add(c);
    }
    
    public Iterator<WeightMatrixHit> execute(Region r) { 
    	LinkedList<WeightMatrixHit> hits = new LinkedList<WeightMatrixHit>();

        String seq = seqgen.execute(r);
        float[] scores = null;
        int s = r.getStart();

        try { 
            for(int k = 0; k < matrices.size(); k++) {
                WeightMatrix matrix = matrices.get(k);
                int width = matrix.matrix.length;
                float cutoff = cutoffs.get(k);

                scores = score(matrix, seq.toCharArray());

                for(int i = 0; i < scores.length; i++) { 
                    if(scores[i] >= cutoff) { 
                        WeightMatrixHit hit = new WeightMatrixHit(r.getGenome(), r.getChrom(),
                                s+i, s+i+width-1, (double)scores[i], '+', matrix);
                        hits.add(hit);
                    }
                }
            }

            seq = SequenceUtils.reverseComplement(seq);

            for(int k = 0; k < matrices.size(); k++) { 
                WeightMatrix matrix = matrices.get(k);
                int width = matrix.matrix.length;
                float cutoff = cutoffs.get(k);
                scores = score(matrix, seq.toCharArray());

                // i = seq.length - ip - width + 1

                for(int ip = 0; ip < scores.length; ip++) { 
                    if(scores[ip] >= cutoff) {
                        int i = scores.length - ip - width;
                        WeightMatrixHit hit = new WeightMatrixHit(r.getGenome(), r.getChrom(),
                                s+i, s+i+width-1, (double)scores[ip], '-', matrix);
                        hits.add(hit);
                    }
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) { 
            e.printStackTrace(System.err);
        }
       	
    	return hits.iterator();
    }

    public float[] score(WeightMatrix matrix, char[] sequence) {
        float[] results = new float[sequence.length];
        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i < sequence.length - length+1; i++) {
            float score = (float)0.0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            results[i] = score;
        }
        for (int i = sequence.length - length+1; i < sequence.length; i++) {
            results[i] = (float)0.0;
        }

        return results;
    }

}
