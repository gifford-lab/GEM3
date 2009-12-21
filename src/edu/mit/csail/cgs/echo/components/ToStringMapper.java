package edu.mit.csail.cgs.echo.components;

import java.util.Map;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.ewok.verbs.Mapper;

public class ToStringMapper implements SelfDescribingVerb, Mapper<Object,String> {
	
	public ToStringMapper() {}
	
	private static final String[] inputNames = { "Objects" };
	private static final EchoType[] inputClasses = { EchoType.OBJECT_CLASS };

	public EchoType[] getInputClasses() { return inputClasses; }
	public String[] getInputNames() { return inputNames; }

	public EchoType getOutputClass() {
		return new ClassType(String.class);
	}

	public EchoType[] getParameterClasses() {
		return null;
	}

	public String[] getParameterNames() {
		return null;
	}

	public void init(Map<String, Object> params) {
	}

	public String execute(Object a) {
		return a.toString();
	}
}
