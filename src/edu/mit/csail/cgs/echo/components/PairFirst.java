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

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Pair;

public class PairFirst<X,Y> 
    implements Mapper<Pair<X,Y>,X>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private static final String[] inputNames = { "Pairs" };
    private static final EchoType[] inputTypes = { new PairType(EchoType.OBJECT_CLASS, EchoType.OBJECT_CLASS) };
    
    private EchoType outputType;
    
    public PairFirst() { 
        outputType = EchoType.OBJECT_CLASS;
    }

    public X execute(Pair<X, Y> a) {
        return a.getFirst();
    }

    public EchoType[] getInputClasses() { return inputTypes; }
    public String[] getInputNames() { return inputNames; }

    public EchoType getOutputClass() { return outputType; }

    public EchoType[] getParameterClasses() { return null; }
    public String[] getParameterNames() { return null; }

    public void init(Map<String, Object> params) {}

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            outputType = EchoType.OBJECT_CLASS;
        }
    }

    public void clearParameter(String n) {}

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            if(c instanceof PairType) { 
                PairType pt = (PairType)c;
                outputType = pt.getFirstType();
            }
        }
    }

    public void setParameter(String n, EchoType c) {}

}
