/*
 * Created on Apr 13, 2007
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import edu.mit.csail.cgs.ewok.types.*;

public class UniqueFilter<X> implements Filter<X,X>, DependentSelfDescribingVerb {

	private Set<X> previousValues;
	private EchoType outputClass;

	public UniqueFilter() {
		previousValues = new HashSet<X>();
		outputClass = EchoType.OBJECT_CLASS;
	}

	public void reset() { 
		previousValues.clear();
	}

	public X execute(X a) {
		if(previousValues.contains(a)) { 
			return null;
		} else { 
			previousValues.add(a);
			return a;
		}
	}

	private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS };
	private static final String[] inputNames = { "Objects" };

	private static final EchoType[] paramClasses = null;
	private static final String[] paramNames = null;

	public EchoType[] getInputClasses() { return inputClasses; }
	public String[] getInputNames() { return inputNames; }
	public EchoType getOutputClass() { return outputClass; }
	public EchoType[] getParameterClasses() { return paramClasses; }
	public String[] getParameterNames() { return paramNames; }

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

	public void init(Map<String, Object> params) {
		reset();
	}
}
