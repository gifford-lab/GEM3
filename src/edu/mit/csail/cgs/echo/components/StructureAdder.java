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

public class StructureAdder<X>
    implements Mapper<X,Map<String,Object>>, SelfDescribingVerb, DependentSelfDescribingVerb {
    
    private static final String[] inputNames = { "Structures", "Objects" };
    private static final EchoType[] constraintTypes = { new StructuredType("structure"), EchoType.OBJECT_CLASS };
    
    private static final String[] paramNames = null;
    private static final EchoType[] paramTypes = null;

    private EchoType outputType; 
    private EchoType[] inputTypes;
    private String key;
    
    public StructureAdder() { 
        key = "unknown";
        inputTypes = new EchoType[2];
        inputTypes[0] = new StructuredType("structure");
        inputTypes[1] = EchoType.OBJECT_CLASS;
        outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
    }

    public StructureAdder(String k) { 
        key = k;
        inputTypes = new EchoType[2];
        inputTypes[0] = new StructuredType("structure");
        inputTypes[1] = EchoType.OBJECT_CLASS;
        outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
    }

    public Map execute(X a) {
        HashMap<String,Object> map = (HashMap<String,Object>)a;
        map.put(key, a);
        return map;
    }

    public EchoType[] getInputClasses() { return constraintTypes; }
    public String[] getInputNames() { return inputNames; }
    public EchoType getOutputClass() { return outputType; }
    public EchoType[] getParameterClasses() { return paramTypes; }
    public String[] getParameterNames() { return paramNames; }

    public void init(Map<String, Object> params) {
    }

    public void clearInput(String n) {
        if(n.equals(inputNames[0])) { 
            inputTypes[0] = new StructuredType("structure");
            outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
        }
        if(n.equals(inputNames[1])) { 
            inputTypes[1] = EchoType.OBJECT_CLASS;
            outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
        }
    }

    public void setInput(String n, EchoType c) {
        if(n.equals(inputNames[0])) {
            if(c instanceof StructuredType) {
                StructuredType st = (StructuredType)c;
                inputTypes[0] = st;
                outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
            }
        }
        if(n.equals(inputNames[1])) { 
            inputTypes[1] = c;
            outputType = new StructuredType((StructuredType)inputTypes[0], key, inputTypes[1]);
        }
    }

    public void clearParameter(String n) {
    }
    
    public void setParameter(String n, EchoType c) {
    }
}
