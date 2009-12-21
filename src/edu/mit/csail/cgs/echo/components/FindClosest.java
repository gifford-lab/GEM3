/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.Collection;
import java.util.Map;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Pair;

public class FindClosest 
    implements Mapper<Pair<Point,Collection<Region>>,Region>, SelfDescribingVerb {
    
    private static final String[] inputNames = { "RegionPairs" };
    private static final EchoType[] inputClasses = { 
        new PairType(new ClassType(Point.class), new CollectionType(new ClassType(Region.class))) };
    private static final EchoType outputClass = new ClassType(Region.class);

    public Region execute(Pair<Point, Collection<Region>> a) {
        int minDist = -1;
        Point base = a.getFirst();
        Region minReg = null;
        Collection<Region> regs = a.getLast();
        
        for(Region reg : regs) {
            if(reg.getChrom().equals(base.getChrom())) { 
                int dist = reg.distance(base);
                if(minDist == -1 || dist < minDist) { 
                    minReg = reg;
                    minDist = dist;
                }
            }
        }

        return minReg;
    }

    public EchoType[] getInputClasses() {
        return inputClasses;
    }

    public String[] getInputNames() {
        return inputNames;
    }

    public EchoType getOutputClass() {
        return outputClass;
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
