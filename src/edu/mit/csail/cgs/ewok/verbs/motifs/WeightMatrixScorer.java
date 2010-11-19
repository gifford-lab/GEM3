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
        String seq = seqgen.execute(r);
        seq = seq.toUpperCase();
        double[] fscores = null, rscores = null;

        try { 
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

    /**
     * Return the highest scoring sequence in the region
     */
    public String getMaxScoreSequence(Region r, double threshold, int extend){
        String seq = seqgen.execute(r);
        seq = seq.toUpperCase();
        String hit=null;

        int length = matrix.length();
        char[] sequence = seq.toCharArray();
        for (int i = 0; i <= sequence.length - length; i++) {
            float score = (float)0.0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            if (score>threshold){
            	threshold = score;
            	int start = Math.max(0, i-extend);
            	int end = Math.min(sequence.length - length, i+length+extend);
            	hit = seq.substring(start, end);
            }
        }
        seq = SequenceUtils.reverseComplement(seq);
        sequence = seq.toCharArray();
        for (int i = 0; i <= sequence.length - length; i++) {
            float score = (float)0.0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            if (score>threshold){
            	int start = Math.max(0, i-extend);
            	int end = Math.min(sequence.length - length, i+length+extend);
            	hit = seq.substring(start, end);
            }
        }
        
        return hit;
    }
}
