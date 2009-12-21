package edu.mit.csail.cgs.utils.models.data;

public interface Transformation<A,B> {
	public B transform(A v);
	public Class<A> fromClass();
	public Class<B> toClass();
}
