/*
 * Created on Apr 18, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;
import edu.mit.csail.cgs.utils.Pair;

public class CountFractionSink implements Sink<Pair<Integer,Integer>>, SelfDescribingVerb {
    
    private long totalCount;
    private long selectedCount;
    private double fraction;
    
    public CountFractionSink() { 
        totalCount = selectedCount = 0;
        fraction = 0.0;
    }
    
    private static final String[] inputNames = { "Counts" };
    private static final EchoType[] inputTypes = { new PairType(new ClassType(Integer.class), new ClassType(Integer.class)) };

    public EchoType[] getInputClasses() {
        return inputTypes;
    }

    public String[] getInputNames() {
        return inputNames;
    }

    public EchoType getOutputClass() {
        return null;
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

    public void consume(Iterator<Pair<Integer, Integer>> itr) {
        while(itr.hasNext()) { 
            consume(itr.next());
        }
    }

    public void consume(Pair<Integer, Integer> val) {
        selectedCount += val.getFirst();
        totalCount += val.getLast();
    }

    public void finish() {
        fraction = totalCount > 0 ? (double)selectedCount / (double)totalCount : 0.0;
        System.out.println("Selected: " + selectedCount);
        System.out.println("Total: " + totalCount);
        System.out.println("Fraction: " + fraction);
    }

    public void init() {
        totalCount = selectedCount = 0;        
    }

}
