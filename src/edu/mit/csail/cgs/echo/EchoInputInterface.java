/*
 * Created on Feb 19, 2007
 */
package edu.mit.csail.cgs.echo;

/**
 * @author Timothy Danford
 */
public interface EchoInputInterface {    
	public void reset();
	public void evaluate() throws EchoException;
	public EchoProcessor getProcessor();
}
