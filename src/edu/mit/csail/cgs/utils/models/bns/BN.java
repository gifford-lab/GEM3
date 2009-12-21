/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;

import edu.mit.csail.cgs.utils.graphs.DirectedAlgorithms;
import edu.mit.csail.cgs.utils.graphs.DirectedGraph;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.probability.FiniteDistribution;

public class BN<X extends Model> {

	public DirectedGraph graph;
	
	private Map<String,BNVar> vars;
	private Map<String,BNCpd> cpds;
	private DataFrame<X> data;
	
	public BN(BN bn) { 
		data = bn.data;
		vars = new TreeMap<String,BNVar>(bn.vars);
		cpds = new TreeMap<String,BNCpd>(bn.cpds);
		graph = new DirectedGraph(bn.graph);
	}
	
	public BN(DataFrame<X> d, String... varNames) { 
		data = d;
		graph = new DirectedGraph();
		vars = new TreeMap<String,BNVar>();
		cpds = new TreeMap<String,BNCpd>();
		
		for(int i = 0; i < varNames.length; i++) { 
			if(vars.containsKey(varNames[i])) { 
				throw new IllegalArgumentException(varNames[i]);
			}
			BNVar var = new BNVar(varNames[i], data.fieldValues(varNames[i]));
			vars.put(varNames[i], var);
			graph.addVertex(varNames[i]);
		}
		
		// No call to learnCPDs(), since we presume that the user will add some edges to 
		// the graph first, and then call that method him/herself.
	}
	
	public BN(DataFrame<X> d, DirectedGraph g) { 
		data = d;
		graph = g;
		vars = new TreeMap<String,BNVar>();
		cpds = new TreeMap<String,BNCpd>();
	
		Set<String> verts = graph.getVertices();

		for(String varName : verts) { 
			BNVar var = new BNVar(varName, data.fieldValues(varName));
			vars.put(varName, var);
		}
		
		learnCPDs();
	}
	
	public FiniteDistribution posterior(X m, String var) { 
		BNVar bnvar = vars.get(var);
		Object bnvalue = bnvar.findValue(m);
		if(bnvalue != null) { 
			int code = bnvar.encode(bnvalue);
			return new FiniteDistribution(bnvar.size(), code);
		}
		
		
		
		return null;
	}
	
	public void print() { 
		DirectedAlgorithms algos = new DirectedAlgorithms(graph);
		Vector<String> verts = algos.getTopologicalOrdering();
		graph.printGraph(System.out);
		for(String name : verts) { 
			cpds.get(name).print();
			System.out.println();
		}
	}
	
	public DataFrame<X> getData() { return data; }
	
	public Set<String> varNames() { return graph.getVertices(); }
	
	public X sample() { 
		try {
			X model = data.getModelClass().newInstance();
			DirectedAlgorithms algos = new DirectedAlgorithms(graph);
			Vector<String> verts = algos.getTopologicalOrdering();

			for(String v : verts) { 
				cpds.get(v).resample(model);
			}
			return model;
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public double logLikelihood(Model m) { 
		double sum = 0.0;
		for(String varName : cpds.keySet()) { 
			sum += cpds.get(varName).logLikelihood(m);
		}
		return sum;
	}
	
	public double logLikelihood(Iterator<? extends Model> ms) { 
		double sum = 0.0;
		while(ms.hasNext()) { 
			sum += logLikelihood(ms.next());
		}
		return sum;
	}
	
	public double logLikelihood() { 
		return logLikelihood(data.iterator());
	}
	
	public BNVar getVar(String v) { 
		return vars.get(v);
	}
	
	public BNCpd getCPD(String v) { 
		return cpds.get(v);
	}
	
	public int countParameters() {
		int count = 0;
		for(String key : cpds.keySet()) { 
			count += cpds.get(key).countParameters();
		}
		return count;
	}
	
	public int countParameters(String child, String... parents) { 
		BNVar cvar = vars.get(child);
		int count = cvar.size();
		
		for(int i = 0; i < parents.length; i++) { 
			count *= vars.get(parents[i]).size();
		}
		
		return count;
	}
	
	public void learnCPDs() { 
		DirectedAlgorithms algos = new DirectedAlgorithms(graph);
		if(algos.hasCycle()) { 
			throw new IllegalStateException("Graph has a cycle.");
		}
		
		for(String name : vars.keySet()) { 
			cpds.put(name, learnCPD(name));
		}		
	}
	
	private BNCpd learnCPD(String node) {
		Set<String> parentNames = graph.getParents(node);
		
		BNVar[] parents = new BNVar[parentNames.size()];
		BNVar child = vars.get(node);
		int pi = 0;
		for(String p : parentNames) { 
			parents[pi++] = vars.get(p);
		}
		
		BNCpd cpd = new BNCpd(parents, child);
		
		cpd.learn(data.iterator());
		
		return cpd;
	}
}
