/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.Map;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class PointToRegion implements SelfDescribingVerb, Mapper<Point,Region> {
    
    public PointToRegion() {}
    
    private static final String[] inputNames = { "Points" };
    private static final EchoType[] inputTypes = { new ClassType(Point.class) };
    private static final EchoType outputType = new ClassType(Region.class);

    public EchoType[] getInputClasses() { return inputTypes; }

    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() {
        return null;
    }

    public String[] getParameterNames() {
        return null;
    }

    public void init(Map<String, Object> params) {
    }

    public Region execute(Point a) {
        return new Region(a.getGenome(), a.getChrom(), a.getLocation(), a.getLocation());
    }
}
