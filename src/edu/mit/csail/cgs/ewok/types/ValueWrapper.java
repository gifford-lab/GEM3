package edu.mit.csail.cgs.ewok.types;

import javax.swing.*;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.NotFoundException;

public class ValueWrapper implements SelfDescribingConstant { 
	private EchoType vclass;
	private Object value;

	public ValueWrapper(Object v) { 
		value = v;
		vclass = new ClassType(v.getClass());
	}

	public Object getConstantValue() { return value; }
	public EchoType getConstantClass() { return vclass; }
	public void setConstantValue(Object v) { value = v; }
    
    public String toString() { return value.toString(); }

}
