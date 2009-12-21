/*
 * Author: tdanford
 * Date: Aug 28, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

/**
 * ATransformation is an "abstract Transformation" -- it's an implementation of the Transformation
 * interface, that automatically fills in the two typing methods.  As such, it's useful for 
 * extending as an anonymous inline class (i.e.
 * <code>
 * Transformation<Foo,Bar> t = new ATransformation<Foo,Bar>(Foo.class, Bar.class) { 
 * 	public Bar transform (Foo a) { 
 * 		...
 * 	}
 * };
 * </code>
 * 
 * @author tdanford
 *
 * @param <X>
 * @param <Y>
 */
public abstract class ATransformation<X,Y> implements Transformation<X,Y> {
	
	private Class<X> from;
	private Class<Y> to;
	
	public ATransformation(Class<X> f, Class<Y> t) { 
		from = f;
		to = t;
	}

	public Class<X> fromClass() { return from; }
	public Class<Y> toClass() { return to; }
}
