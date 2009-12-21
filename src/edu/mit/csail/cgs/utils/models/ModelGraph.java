/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.GenSym;
import edu.mit.csail.cgs.utils.graphs.*;

public class ModelGraph {

	public DirectedGraph graph;
	private Map<Model,String> modelToName;
	private Map<String,Model> nameToModel;
	private GenSym gensym;
	private int maxDepth;
	
	public ModelGraph() { 
		gensym = new GenSym("m");
		graph = new DirectedGraph();
		modelToName = new HashMap<Model,String>();
		nameToModel = new TreeMap<String,Model>();
		maxDepth = -1;
	}
	
	public void setMaxDepth(Integer md) { maxDepth = md; }
	public Integer getMaxDepth() { return maxDepth; }
	
	public Model findModel(String graphName) { return nameToModel.get(graphName); }
	
	public void addModels(Iterator<Model> ms) { 
		while(ms.hasNext()) { 
			addModel(ms.next());
		}
	}
	
	public String addModel(Model m) { 
		return addModel(m, 0);
	}
	
	private String addModel(Model m, int depth) { 
		if(maxDepth != -1 && depth > maxDepth) { 
			return null;
		}
		
		if(!modelToName.containsKey(m)) { 
			String name = gensym.generateSymbol();
			modelToName.put(m, name);
			nameToModel.put(name, m);
			
			graph.addVertex(name);
			
			ModelFieldAnalysis analysis = new ModelFieldAnalysis(m.getClass());
			Vector<Field> fs = analysis.getFields();
			for(Field f : fs) {
				if(Model.isSubclass(f.getType(), Model.class)) { 
					try {
						Model value = (Model)f.get(m);
						if(value != null) {
							String valueName = addModel(value, depth+1);
							if(valueName != null) { 
								graph.addEdge(name, valueName);
							}
						}
					} catch (IllegalAccessException e) {
						e.printStackTrace();
					}
				}
			}
		
			return name;
		} else { 
			return modelToName.get(m);
		}
	}
}
