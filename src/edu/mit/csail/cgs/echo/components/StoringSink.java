/*
 * Created on Feb 22, 2007
 */
package edu.mit.csail.cgs.echo.components;

import java.util.*;
import edu.mit.csail.cgs.ewok.verbs.Sink;
import edu.mit.csail.cgs.ewok.types.*;

public class StoringSink<X> implements SelfDescribingOutput<X>, Sink<X>, DependentSelfDescribingVerb {
    
    private LinkedList<X> values;
	private EchoType outputClass;
    
    public StoringSink() { 
		outputClass = EchoType.OBJECT_CLASS;
		init();
    }

    private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS };
    private static final String[] inputNames = { "Objects" };
    
    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }
    public EchoType getOutputClass() { return outputClass; }
    public EchoType[] getParameterClasses() { return null; }
    public String[] getParameterNames() { return null; }

	public void setInput(String n, EchoType c) { 
		if(n.equals(inputNames[0])) {
			outputClass = c;
		}
	}

	public void clearInput(String n) { 
		if(n.equals(inputNames[0])) { 
			outputClass = EchoType.OBJECT_CLASS;
		}
	}

	public void setParameter(String n, EchoType c) {}
	public void clearParameter(String n) {}

    public void init(Map<String,Object> params) {
		init();
    }

	public void init() {
        values = new LinkedList<X>();
	}
    
    public Collection<X> getValues() { return values; }

    public void consume(Iterator<X> itr) {
        while(itr.hasNext()) { 
			consume(itr.next());
        }
    }

	public void consume(X val) { 
		values.addLast(val);
	}

	public void finish() {}
}
