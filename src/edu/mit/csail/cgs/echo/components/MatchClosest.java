/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.*;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.MultiSink;
import edu.mit.csail.cgs.utils.Pair;

public class MatchClosest implements MultiSink, SelfDescribingVerb, SelfDescribingOutput<Pair<Region,Region>> {
    
    private Set<Region> points;
    private HashMap<String,Set<Region>> regions;
    private LinkedList<Pair<Region,Region>> matching;
    
    public MatchClosest() { 
        points = new HashSet<Region>();
        regions = new HashMap<String,Set<Region>>();
        matching = new LinkedList<Pair<Region,Region>>();
    }
    
    private static final String[] inputNames = { "Sources", "Targets" };
    private static final EchoType[] inputTypes = 
        { new ClassType(Region.class), new ClassType(Region.class) };
    private static final EchoType outputType = new PairType(new ClassType(Region.class), new ClassType(Region.class)); 

    public Collection<Pair<Region, Region>> getValues() {
        return new LinkedList<Pair<Region,Region>>(matching);
    }

    public void consume(String n, Object val) {
        if(n.equals(inputNames[0])) { points.add((Region)val); }
        
        if(n.equals(inputNames[1])) {
            Region r = (Region)val;
            if(!regions.containsKey(r.getChrom())) { 
                regions.put(r.getChrom(), new HashSet<Region>());
            }
            regions.get(r.getChrom()).add(r);
        }
    }

    public void consume(String n, Iterator itr) {
        while(itr.hasNext()) { consume(n, itr.next()); }
    }

    public void finish() {
        doMatching();
    }
    
    public void doMatching() { 
        for(Region p : points) { 
            Region r = findClosestRegion(p);
            Pair<Region,Region> pair = new Pair<Region,Region>(p, r);
            matching.addLast(pair);
        }
    }
    
    public Region findClosestRegion(Region p) { 
        Region closest = null;
        int minDist = -1;

        if(regions.containsKey(p.getChrom())) { 
            for(Region r : regions.get(p.getChrom())) { 
                if(r.getGenome().equals(p.getGenome())) { 
                    int dist = r.distance(p);
                    if(minDist == -1 || dist < minDist) { 
                        closest = r;
                        minDist = dist;
                    }
                }
            }
        }
        
        return closest;
    }
    
    public void init() {
        points.clear();
        regions.clear();
        matching.clear();
    }
    

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
        init();
    }

}
