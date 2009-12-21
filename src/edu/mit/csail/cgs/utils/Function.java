package edu.mit.csail.cgs.utils;

public interface Function<X,Y> {
	public boolean inDomain(X value);
	public Y valueAt(X value);
	
	public static class Constant<A,B> implements Function<A,B> { 
		private B value;
		public Constant(B v) { value = v; }
		
		public boolean inDomain(A value) {
			return true;
		}
		
		public B valueAt(A v) {
			return value;
		}
	}
}
