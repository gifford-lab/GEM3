package edu.mit.csail.cgs.utils.graphs;

import java.util.*;

public class VertexSet {

	private TreeSet<String> verts;
	
	public VertexSet(String vs) { 
		String[] array = vs.split(",");
		verts = new TreeSet<String>();
		for(int i = 0; i < array.length; i++) { 
			verts.add(array[i]);
		}
	}
	
	public VertexSet(Collection<String> v) { 
		verts = new TreeSet<String>(v);
	}
	
	public int size() { return verts.size(); }
	
	public Collection<String> getVertices() { return new LinkedList<String>(verts); }
	
	public int hashCode() { 
		int code = 17;
		for(String v : verts) { 
			code += v.hashCode(); code *= 37;
		}
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof VertexSet)) { return false; }
		VertexSet vs = (VertexSet)o;
		if(verts.size() != vs.verts.size()) { 
			return false;
		}
		
		for(String v : verts) { if(!vs.verts.contains(v)) { return false; } }
		
		return true;
	}
	
	public boolean isClique(Graph g) { 
		for(String v1 : verts) { 
			for(String v2 : verts) { 
				if(!v1.equals(v2)) { 
					if(!g.isNeighbor(v1, v2) ||
							!g.isNeighbor(v2, v1)) { 
						return false;
					}
				}
			}
		}
		return true;
	}
}
