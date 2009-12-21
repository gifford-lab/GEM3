package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class WeightMatrixScorer implements Mapper<Region,WeightMatrixScoreProfile> {

    private WeightMatrix matrix;
    private SequenceGenerator seqgen;
    
    public WeightMatrixScorer(WeightMatrix m) {
    	matrix = m;
        seqgen = new SequenceGenerator();
    }
    
    public WeightMatrixScoreProfile execute(Region r) { 
    	LinkedList<WeightMatrixHit> hits = new LinkedList<WeightMatrixHit>();

        String seq = seqgen.execute(r);
        seq = seq.toUpperCase();
        double[] fscores = null, rscores = null;
        int s = r.getStart();

        try { 
        	int width = matrix.matrix.length;
        	fscores = score(matrix, seq.toCharArray(), '+');
            seq = SequenceUtils.reverseComplement(seq);
            rscores = score(matrix, seq.toCharArray(), '-');

        } catch (ArrayIndexOutOfBoundsException e) { 
            e.printStackTrace(System.err);
        }
       	
        return new WeightMatrixScoreProfile(matrix, fscores, rscores);
    }
    public WeightMatrixScoreProfile execute(String seq) { 
    	seq = seq.toUpperCase();
        double[] fscores = null, rscores = null;
    	fscores = score(matrix, seq.toCharArray(), '+');
        seq = SequenceUtils.reverseComplement(seq);
        rscores = score(matrix, seq.toCharArray(), '-');
        return new WeightMatrixScoreProfile(matrix, fscores, rscores);
    }

    public double[] score(WeightMatrix matrix, char[] sequence, char strand) {
        double[] results = new double[sequence.length];
        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i <= sequence.length - length; i++) {
            float score = (float)0.0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            
            if(strand=='-') { 
            	results[sequence.length-length-i] = score;
            } else { 
            	results[i] = score;
            }
        }
       	for (int i = (sequence.length - length)+1; i < sequence.length; i++) {
       		results[i] = (float)matrix.getMinScore();
       	}
       
        return results;
    }

}
