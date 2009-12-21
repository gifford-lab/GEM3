package edu.mit.csail.cgs.utils.graphs.layout;

import edu.mit.csail.cgs.utils.graphs.*;
import java.io.*;
import java.util.*;
import java.awt.Point;

public interface LayoutEngine<G extends Graph> { 
	public GraphLayout layoutGraph(G graph);
	public Set<String> getParameters();
	public void setParameter(String key, String value);
}
