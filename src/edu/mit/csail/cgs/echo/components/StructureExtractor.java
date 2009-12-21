/*
 * Created on Apr 18, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.*;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class StructureExtractor<X>
    implements Mapper<Map<String,Object>,X>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private static final String[] inputNames = { "Structures" };
    
    private static final String[] paramNames = null;
    private static final EchoType[] paramTypes = null;

    private EchoType outputType;
    private EchoType[] inputTypes;
    private String key;
    
    public StructureExtractor() { 
        key = "unknown";
        inputTypes = new EchoType[1];
        inputTypes[0] = new StructuredType("structure", key, EchoType.OBJECT_CLASS);
        outputType = EchoType.OBJECT_CLASS;
    }

    public StructureExtractor(String k) { 
        key = k;
        inputTypes = new EchoType[1];
        inputTypes[0] = new StructuredType("structure", key, EchoType.OBJECT_CLASS);
        outputType = EchoType.OBJECT_CLASS;
    }

    public X execute(Map<String,Object> map) {
        return (X)map.get(key);
    }

    public EchoType[] getInputClasses() { return inputTypes; }
    public String[] getInputNames() { return inputNames; }
    public EchoType getOutputClass() { return outputType; }
    public EchoType[] getParameterClasses() { return paramTypes; }
    public String[] getParameterNames() { return paramNames; }

    public void init(Map<String, Object> params) {
    }

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) {
            outputType = EchoType.OBJECT_CLASS;
        }
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) {
            if(c instanceof StructuredType) { 
                StructuredType st = (StructuredType)c;
                outputType = st.getTypeValue(key);
            }
        }
    }

    public void clearParameter(String n) {
    }
    
    public void setParameter(String n, EchoType c) {
    }
}
