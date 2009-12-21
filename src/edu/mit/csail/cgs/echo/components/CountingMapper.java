package edu.mit.csail.cgs.echo.components;

import java.util.Iterator;
import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class CountingMapper<X> 
	implements Mapper<Iterator<X>,Integer>, SelfDescribingVerb {
	
	public CountingMapper() {}

	public Integer execute(Iterator<X> a) {
		int count = 0;
		while(a.hasNext()) { 
			a.next();
			count += 1;
		}
		return count;
	}
	
	private EchoType[] inputClasses = { new SequenceType(EchoType.OBJECT_CLASS) };
	private String[] inputNames = { "Objects" };

	public EchoType[] getInputClasses() { return inputClasses; }
	public String[] getInputNames() { return inputNames; }

	public EchoType getOutputClass() {
		return new ClassType(Integer.class);
	}

	public EchoType[] getParameterClasses() {
		return null;
	}

	public String[] getParameterNames() {
		return null;
	}

	public void init(Map<String, Object> params) {
	}
}
