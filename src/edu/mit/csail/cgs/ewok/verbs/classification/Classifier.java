/*
 * Created on Mar 15, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.classification;

import edu.mit.csail.cgs.ewok.verbs.Mapper;

/**
 * @author tdanford
 */
public interface Classifier<X> extends Mapper<X,Integer> {
    public int getNumClasses();
}
