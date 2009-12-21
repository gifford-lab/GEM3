package edu.mit.csail.cgs.ewok.verbs;

public interface BiCombiner<X,Y,Z> {
	public Z execute(X a, Y b);
}
