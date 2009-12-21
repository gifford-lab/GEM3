/*
 * Created on Feb 16, 2007
 */
package edu.mit.csail.cgs.ewok.types;



public interface SelfDescribingVerb extends SelfDescribingParameterized {    
    public String[] getInputNames();
    public EchoType[] getInputClasses();
    public EchoType getOutputClass();
}
