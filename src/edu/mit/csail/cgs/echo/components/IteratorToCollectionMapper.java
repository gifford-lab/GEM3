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
import java.util.Iterator;
import java.util.Map;
import java.util.LinkedList;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class IteratorToCollectionMapper<X> implements 
    Mapper<Iterator<X>,Collection<X>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private EchoType outputType;
    
    private static final EchoType[] inputTypes = { new SequenceType(EchoType.OBJECT_CLASS) };
    private static final String[] inputNames = { "Iterator" };
    
    public IteratorToCollectionMapper() { 
        outputType = new CollectionType(EchoType.OBJECT_CLASS);
    }

    public Collection<X> execute(Iterator<X> a) {
        LinkedList<X> vals = new LinkedList<X>();
        while(a.hasNext()) { vals.addLast(a.next()); }
        return vals;
    }

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

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            outputType = new CollectionType(EchoType.OBJECT_CLASS);
        }
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            if(c instanceof SequenceType) { 
                SequenceType st = (SequenceType)c;
                outputType = new CollectionType(st.getInnerType());
            } else { 
                outputType = new CollectionType(c);
            }
        }
    }

    public void clearParameter(String n) {
    }

    public void setParameter(String n, EchoType c) {
    }

}
