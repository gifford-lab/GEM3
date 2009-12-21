package edu.mit.csail.cgs.utils.graphs;

public interface CycleChecker {
	public boolean containsCycle();
	public boolean checkCycleWithEdge(String v1, String v2);
}