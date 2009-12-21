/*
 * Author: tdanford
 * Date: Sep 29, 2008
 */
package edu.mit.csail.cgs.utils.models;


public interface Modeleable {
	public Class getModelClass();
	public Model asModel();
}
