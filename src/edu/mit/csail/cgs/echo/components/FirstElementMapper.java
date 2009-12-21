/*
 * Created on Apr 17, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.Map;

import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.echo.*;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.Closeable;

public class FirstElementMapper<X> implements SelfDescribingVerb, DependentSelfDescribingVerb, Mapper<Iterator<X>,X> {
    
    private EchoType outputType;
    private static final EchoType[] inputTypes = { new SequenceType(EchoType.OBJECT_CLASS) };
    private static final String[] inputNames = { "Iterator" };
    
    public FirstElementMapper() {
        outputType = EchoType.OBJECT_CLASS;
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

    public void clearInput(String n) {
        outputType = EchoType.OBJECT_CLASS;
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            if(c instanceof SequenceType) { 
                SequenceType st = (SequenceType)c;
                outputType = st.getInnerType();
            } else { 
                outputType = c;
            }
        }
    }

    public void clearParameter(String n) {}

    public void setParameter(String n, EchoType c) {}

    public EchoType[] getParameterClasses() { return null; }
    public String[] getParameterNames() { return null; }

    public void init(Map<String, Object> params) {}

    public X execute(Iterator<X> a) {
        X val = a.hasNext() ? a.next() : null;
        if(a instanceof Closeable) { ((Closeable)a).close(); }
        return val;
    }

}
