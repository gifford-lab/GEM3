package edu.mit.csail.cgs.echo.gui;

import java.util.*;

public interface EchoEdgeConnector {
	public Collection<EchoEdge> getEdges();
	public void connect(EchoComponent comp);
	public void clearEdges();
	public void clearEdge(EchoComponent c);
}
