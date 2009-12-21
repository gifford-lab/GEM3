/*
 * Author: tdanford
 * Date: Jan 26, 2009
 */
package edu.mit.csail.cgs.utils.models;

public interface ModelListener<T extends Model> {
	public void modelChanged(T model);
}
