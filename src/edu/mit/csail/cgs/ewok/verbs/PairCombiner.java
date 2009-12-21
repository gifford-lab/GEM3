package edu.mit.csail.cgs.ewok.verbs;

import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.Pair;

public class PairCombiner<X,Y> implements BiCombiner<X,Y,Pair<X,Y>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private EchoType outputClass;
    private EchoType[] currentInputs;
	
	public PairCombiner() {
        currentInputs = new EchoType[2];
        currentInputs[0] = currentInputs[1] = EchoType.OBJECT_CLASS;
        outputClass = new PairType(EchoType.OBJECT_CLASS, EchoType.OBJECT_CLASS);
    }
    
    private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS, EchoType.OBJECT_CLASS };
    private static final String[] inputNames = { "First", "Second" };

	public Pair<X,Y> execute(X a, Y b) { 
		return new Pair<X,Y>(a, b);
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

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            currentInputs[0] = EchoType.OBJECT_CLASS;
        }
        if(n.equals(inputNames[1])) { 
            currentInputs[1] = EchoType.OBJECT_CLASS;
        }
        outputClass = new PairType(currentInputs[0], currentInputs[1]);
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            currentInputs[0] = c;
        }
        if(n.equals(inputNames[1])) { 
            currentInputs[1] = c;
        }
        outputClass = new PairType(currentInputs[0], currentInputs[1]);
    }

    public void clearParameter(String n) {
    }

    public void setParameter(String n, EchoType c) {
    }
}
