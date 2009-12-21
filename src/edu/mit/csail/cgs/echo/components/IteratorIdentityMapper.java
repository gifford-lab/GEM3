/*
 * Created on Apr 17, 2007
 */
package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class IteratorIdentityMapper<X> 
    implements Mapper<Iterator<X>,Iterator<X>>, SelfDescribingVerb, DependentSelfDescribingVerb {

    private static final String[] inputNames = { "Iterator" };
    private static final EchoType[] inputClasses = { new SequenceType(EchoType.OBJECT_CLASS) };
    
    private EchoType outputClass;

    public IteratorIdentityMapper() {
        outputClass = new SequenceType(EchoType.OBJECT_CLASS);
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

    public Iterator<X> execute(Iterator<X> a) {
        return a;
    }

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            outputClass = new SequenceType(EchoType.OBJECT_CLASS);
        }
    }

    public void clearParameter(String n) {}

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) {
            if(c instanceof SequenceType) { 
                outputClass = c;
            } else { 
                outputClass = new SequenceType(c);
            }
        }
    }

    public void setParameter(String n, EchoType c) {
    } 
}
