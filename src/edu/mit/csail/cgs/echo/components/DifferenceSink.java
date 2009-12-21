/*
 * Created on Feb 22, 2007
 */
package edu.mit.csail.cgs.echo.components;

import java.util.*;
import edu.mit.csail.cgs.echo.Reverb;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.ewok.verbs.MultiSink;

public class DifferenceSink<X> implements SelfDescribingOutput<X>, DependentSelfDescribingVerb, MultiSink<X> {
    
    private Set<X> firstValues, secondValues;
	private EchoType outputClass;
    private EchoType[] inputClasses;
    
    public DifferenceSink() { 
		outputClass = EchoType.OBJECT_CLASS;
		inputClasses = new EchoType[2];
		inputClasses[0] = EchoType.OBJECT_CLASS;
		inputClasses[1] = EchoType.OBJECT_CLASS;
		init();
    }

	public static Class getLUBClass(Class c1, Class c2) { 
		if(Reverb.isSubclass(c1, c2)) { return c2; }
		if(Reverb.isSubclass(c2, c1)) { return c1; }
		return Object.class;
	}

	public static EchoType getLUBType(EchoType c1, EchoType c2) { 
		if(c1.isSubType(c2)) { return c2; }
		if(c2.isSubType(c1)) { return c1; }
		return EchoType.OBJECT_CLASS;
	}

    private static final String[] inputNames = { "First", "Second" };
    
    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }
    public EchoType getOutputClass() { return outputClass; }
    public EchoType[] getParameterClasses() { return null; }
    public String[] getParameterNames() { return null; }

	public void setInput(String n, EchoType c) { 
		if(n.equals(inputNames[0])) { 
			EchoType lub = getLUBType(c, inputClasses[1]);
			outputClass = inputClasses[0] = inputClasses[1] = lub;
		}
		if(n.equals(inputNames[1])) { 
			EchoType lub = getLUBType(c, inputClasses[0]);
			outputClass = inputClasses[0] = inputClasses[1] = lub;
		}
	}

	public void clearInput(String n) { setInput(n, EchoType.OBJECT_CLASS); }

	public void setParameter(String n, EchoType c) {}
	public void clearParameter(String n) {}

    public void init(Map<String,Object> params) {
		init();
    }

	public void init() {
		firstValues = new HashSet<X>();
		secondValues = new HashSet<X>();
	}
    
    public Collection<X> getValues() { 
		SetTools<X> tools = new SetTools<X>();
		return tools.subtract(firstValues, secondValues);
	}

    public void consume(String n, Iterator<X> itr) {
        while(itr.hasNext()) { 
			consume(n, itr.next());
        }
    }

	public void consume(String n, X val) { 
		if(n.equals(inputNames[0])) { firstValues.add(val); }
		if(n.equals(inputNames[1])) { secondValues.add(val); }
	}

	public void finish() {
	}
}
