package edu.mit.csail.cgs.echo;

public interface InputQueue<X> { 
	public boolean isFinished(); 
	public boolean isEmpty(); 
	public X getFirstValue();
}
