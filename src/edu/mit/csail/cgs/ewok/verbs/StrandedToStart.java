package edu.mit.csail.cgs.ewok.verbs;

import java.util.Map;

import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.ewok.types.*;

/**
 * Maps a StrandedRegion to its start site
 */
public class StrandedToStart<X extends StrandedRegion> implements Mapper<X,Point>, SelfDescribingVerb {

    public StrandedToStart() {}

    public Point execute(X a) {
        switch(a.getStrand()) { 
        case '+':
            return new Point(a.getGenome(),
                             a.getChrom(),
                             a.getStart());
        case '-':
            return new Point(a.getGenome(),
                             a.getChrom(),
                             a.getEnd());
        default:
            throw new IllegalArgumentException("Don't understand strand " + a.getStrand());
        }
    }
    
    
    private static final String[] inputNames = { "StrandedRegions" };
    private static final EchoType[] inputTypes = { new ClassType(StrandedRegion.class) };
    private static final EchoType outputType = new ClassType(Point.class); 
    
    public EchoType[] getInputClasses() {
        return inputTypes;
    }

    public String[] getInputNames() {
        return inputNames;
    }

    public EchoType getOutputClass() {
        return outputType;
    }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }
}
