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
import edu.mit.csail.cgs.ewok.verbs.Filter;
import edu.mit.csail.cgs.utils.Pair;

public class NonNullPairFilter<X,Y> 
    implements Filter<Pair<X,Y>,Pair<X,Y>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private EchoType[] inputTypes;
    private EchoType outputType;
    
    private static final String[] inputNames = { "Pairs" };
    
    public NonNullPairFilter() { 
        inputTypes = new EchoType[1];
        inputTypes[0] = new PairType(EchoType.OBJECT_CLASS, EchoType.OBJECT_CLASS);
        outputType = new PairType(EchoType.OBJECT_CLASS,EchoType.OBJECT_CLASS);
    }

    public Pair<X, Y> execute(Pair<X, Y> a) {
        if(a.getLast() != null) { 
            return a; 
        } else { 
            return null;
        }
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
    }

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            outputType = new PairType(EchoType.OBJECT_CLASS,EchoType.OBJECT_CLASS);            
        }
    }

    public void clearParameter(String n) {
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            if(c instanceof PairType) { 
                outputType = c;
            }
        }
    }

    public void setParameter(String n, EchoType c) {
    } 

}
