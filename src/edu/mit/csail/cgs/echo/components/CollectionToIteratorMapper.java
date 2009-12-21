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

import edu.mit.csail.cgs.ewok.types.CollectionType;
import edu.mit.csail.cgs.ewok.types.DependentSelfDescribingVerb;
import edu.mit.csail.cgs.ewok.types.EchoType;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.ewok.types.SequenceType;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class CollectionToIteratorMapper<X> implements 
    Mapper<Collection<X>,Iterator<X>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private EchoType outputType;
    
    private static final EchoType[] inputTypes = { new CollectionType(EchoType.OBJECT_CLASS) };
    private static final String[] inputNames = { "Collection" };
    
    public CollectionToIteratorMapper() { 
        outputType = new SequenceType(EchoType.OBJECT_CLASS);
    }

    public Iterator<X> execute(Collection<X> a) {
        return a.iterator();
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
            outputType = new SequenceType(EchoType.OBJECT_CLASS);
        }
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            if(c instanceof CollectionType) { 
                CollectionType st = (CollectionType)c;
                outputType = new SequenceType(st.getInnerType());
            }
        }
    }

    public void clearParameter(String n) {
    }

    public void setParameter(String n, EchoType c) {
    }

}
