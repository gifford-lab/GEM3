package edu.mit.csail.cgs.ewok.verbs;

import edu.mit.csail.cgs.datasets.general.Region;

public class DistanceFilter<X extends Region> implements Filter<X,X> {

    private int maxdistance;
    private Region base;

    public DistanceFilter(int d, Region b) {
        maxdistance = d;
        base = b;
    }
    public X execute(X other) {
        if (base.distance(other) < maxdistance) {
            return other;
        } else {
            return null;
        }
    }
}