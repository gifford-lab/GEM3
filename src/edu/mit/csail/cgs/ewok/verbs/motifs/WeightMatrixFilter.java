package edu.mit.csail.cgs.ewok.verbs.motifs;

import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Filter;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class WeightMatrixFilter implements Filter<Region,Region> {

    private Vector<WeightMatrix> matrices;
    private Vector<Float> cutoffs;
    private SequenceGenerator seqgen;
    
    public WeightMatrixFilter(WeightMatrix m, float c) {
        matrices = new Vector<WeightMatrix>();
        cutoffs = new Vector<Float>();
        matrices.add(m);
        cutoffs.add(c);
        seqgen = new SequenceGenerator();
    }
    
    public WeightMatrixFilter() {
        matrices = new Vector<WeightMatrix>();
        cutoffs = new Vector<Float>();
        seqgen = new SequenceGenerator();
    }
    
    public void addWeightMatrix(WeightMatrix m, float c) { 
        matrices.add(m);
        cutoffs.add(c);
    }
    
    public Region execute(Region r) { 
        String seq = seqgen.execute(r);
        float[] scores = null;
        int s = r.getStart();

        try { 
            for(int k = 0; k < matrices.size(); k++) {
                WeightMatrix matrix = matrices.get(k);
                int width = matrix.matrix.length;
                float cutoff = cutoffs.get(k);

                if(score(matrix, seq.toCharArray(), cutoff)) { 
                	return r;
                }
            }

            seq = SequenceUtils.reverseComplement(seq);

            for(int k = 0; k < matrices.size(); k++) { 
                WeightMatrix matrix = matrices.get(k);
                int width = matrix.matrix.length;
                float cutoff = cutoffs.get(k);
                if(score(matrix, seq.toCharArray(), cutoff)) { 
                	return r;
                }
            }
        } catch (ArrayIndexOutOfBoundsException e) { 
            e.printStackTrace(System.err);
        }
       	
    	return null;
    }

    public boolean score(WeightMatrix matrix, char[] sequence, float thresh) {
        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i < sequence.length - length; i++) {
            float score = (float)0.0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            if(score >= thresh) { return true; }
        }
        return false;
    }

}
