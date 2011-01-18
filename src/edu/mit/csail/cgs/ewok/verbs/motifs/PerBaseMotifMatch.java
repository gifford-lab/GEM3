package edu.mit.csail.cgs.ewok.verbs.motifs;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.datasets.motifs.*;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class PerBaseMotifMatch implements Mapper<String,Double[]> {
	 
    private WeightMatrix matrix;
    private char[] hitStrands;
    
    public PerBaseMotifMatch(WeightMatrix matrix) {
        this.matrix = matrix;
    }

    public Double[] execute(String sequence) {
        return execute(sequence.toCharArray());
    }
    public Double[] execute(char[] sequence) {
        Double[] results = new Double[sequence.length];
        hitStrands = new char[sequence.length];

        /* scan through the sequence */
        int length = matrix.length();
        for (int i = 0; i < sequence.length - length; i++) {
            double score = 0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            results[i] = score;
            hitStrands[i]='+';
        }
        for (int i = sequence.length - length; i < sequence.length; i++) {
            results[i] = 0.0;
        }
        /* now reverse complement the sequence and scan that.
           basically the same as before but with slightly different
           computation of the hit coordinates */
        SequenceUtils.reverseComplement(sequence);
        for (int i = 0; i < sequence.length - length; i++) {
            double score = 0;
            for (int j = 0; j < length; j++) {
                score += matrix.matrix[j][sequence[i+j]];
            }
            int realindex = sequence.length - i - length;
            if (score > results[realindex]) {
                results[realindex] = score;
                hitStrands[realindex]='-';
            }
        }        
        SequenceUtils.reverseComplement(sequence);
        return results;
    }
    public char [] getHitStrands(){return hitStrands;}
}
