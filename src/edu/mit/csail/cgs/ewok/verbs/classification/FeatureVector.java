/*
 * Created on Mar 15, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.classification;

import edu.mit.csail.cgs.ewok.verbs.Mapper;
import java.util.*;

/**
 * @author tdanford
 */
public class FeatureVector<X> implements Mapper<X,Vector<Double>> {
    
    private Vector<Feature<X>> features;

    public FeatureVector(Vector<Feature<X>> fv) {
        features = new Vector<Feature<X>>(fv);
    }
    
    public int getDimension() { return features.size(); }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public Vector<Double> execute(X value) {
        Vector<Double> vec = new Vector<Double>();
        for(int i = 0; i < features.size(); i++) { 
            vec.add(features.get(i).execute(value));
        }
        return vec;
    }

}
