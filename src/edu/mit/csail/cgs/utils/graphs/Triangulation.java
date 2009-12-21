package edu.mit.csail.cgs.utils.graphs;

public interface Triangulation {
	public UndirectedGraph triangulate(UndirectedGraph ug);
	
	public static class Default implements Triangulation {

		public UndirectedGraph triangulate(UndirectedGraph ug) {
			UndirectedGraph nug = new UndirectedGraph(ug);
			UndirectedCycleChecker cc = new UndirectedCycleChecker(nug);
			
			return nug;
		} 
		
	}
}
