/*
 * Author: tdanford
 * Date: Mar 21, 2009
 */
package edu.mit.csail.cgs.utils.models;

public interface Accessor<T extends Model> {
	public Object get(T object);
	public void set(T object, Object val);
	public Class getType();
	public String getName();
	public String getBaseName();
}
