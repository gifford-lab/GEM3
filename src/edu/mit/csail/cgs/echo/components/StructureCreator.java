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

public class StructureCreator<X>
    implements Mapper<X,Map<String,Object>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private static final String[] inputNames = { "Objects" };
    private static final EchoType[] inputTypes = { EchoType.OBJECT_CLASS };
    
    private static final String[] paramNames = null;
    private static final EchoType[] paramTypes = null;

    private EchoType outputType, inputType;
    private String key;
    
    public StructureCreator() { 
        key = "unknown";
        inputType = EchoType.OBJECT_CLASS;
        outputType = new StructuredType("structure", key, inputType);
    }

    public StructureCreator(String k) { 
        key = k;
        inputType = EchoType.OBJECT_CLASS;
        outputType = new StructuredType("structure", key, inputType);
    }

    public Map execute(X a) {
        HashMap<String,Object> map = new HashMap<String,Object>();
        map.put(key, a);
        return map;
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
            inputType = EchoType.OBJECT_CLASS;
            outputType = new StructuredType("structure", key, inputType);
        }
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) { 
            inputType = c;
            outputType = new StructuredType("structure", key, inputType);
        }
    }

    public void clearParameter(String n) {
    }
    
    public void setParameter(String n, EchoType c) {
    }
}
