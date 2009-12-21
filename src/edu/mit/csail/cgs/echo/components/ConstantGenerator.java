package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Generator;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;
import edu.mit.csail.cgs.utils.iterators.SingleIterator;

public class ConstantGenerator implements Generator, DependentSelfDescribingVerb {
	
	private SelfDescribingConstant konst;
	private EchoType outputClass;
	
	public ConstantGenerator() { 
		konst = null;
		outputClass = null;
	}
	
	public ConstantGenerator(SelfDescribingConstant k) { 
		konst = k;
	}

	public void setInput(String n, EchoType c) {
	}

	public void clearInput(String n) {
	}

	public void setParameter(String n, EchoType c) { 
		if(n.equals(paramNames[0])) { 
			outputClass = c;
		}
	}

	public void clearParameter(String n) { 
		if(n.equals(paramNames[0])) { 
			outputClass = EchoType.OBJECT_CLASS;
		}
	}

	public static final EchoType[] paramClasses = { EchoType.OBJECT_CLASS };
	public static final String[] paramNames = { "Value" };

	public EchoType[] getInputClasses() { return null; }
	public String[] getInputNames() { return null; }

	public EchoType getOutputClass() { return outputClass; }

	public EchoType[] getParameterClasses() { return paramClasses; }
	public String[] getParameterNames() { return paramNames; }

	public void init(Map<String, Object> params) {
		konst = new ValueWrapper(params.get(paramNames[0]));
	}

	public Iterator execute() {
		return konst != null ? 
				new SingleIterator(konst.getConstantValue()) : 
				new EmptyIterator(); 
	}

	public String toString() { return "ConstantGenerator"; }
}
