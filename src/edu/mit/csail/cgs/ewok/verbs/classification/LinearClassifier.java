/*
 * Created on Mar 15, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.classification;

import java.util.*;

/**
 * @author tdanford
 */
public class LinearClassifier implements Classifier<Vector<Double>> {
    
    private Vector<Double> coefficients;

    public LinearClassifier(Vector<Double> coeffs) {
        coefficients = new Vector<Double>(coeffs);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.classification.Classifier#getNumClasses()
     */
    public int getNumClasses() {
        return 2;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public Integer execute(Vector<Double> vec) {
        double sum = 0.0;
        if(vec.size() != coefficients.size()-1) { throw new IllegalArgumentException(); }
        
        for(int i = 0; i < vec.size(); i++) {
            sum += coefficients.get(i) * vec.get(i);
        }
        sum += coefficients.get(coefficients.size()-1);
        
        if(sum >= 0.0) { 
            return 1;
        } else { 
            return 0;
        }
    }

}
