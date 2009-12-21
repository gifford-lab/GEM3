/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.*;

import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.probability.FiniteDistribution;

public class BNCpd {
	
	private BNVar[] parents;
	private BNVar child;

	private Map<BNValues,FiniteDistribution> cpd;
	
	public BNCpd(BNVar[] pars, BNVar chld) { 
		parents = pars.clone();
		child = chld;
		cpd = new LinkedHashMap<BNValues,FiniteDistribution>();
		int childSize = child.size();
		
		BNValuesIterator itr = new BNValuesIterator(parents);
		while(itr.hasNext()) { 
			BNValues parvals = itr.next();
			cpd.put(parvals, new FiniteDistribution(childSize));
		}
	}
	
	public void print() { 
		System.out.println(String.format("CPD: %s", child.getName()));
		for(BNValues vals : cpd.keySet()) { 
			System.out.println(String.format("%s -> %s", vals.toString(), cpd.get(vals).toString()));
		}
	}
	
	public double logLikelihood(Iterator<Model> obs) { 
		double sum = 0.0;
		while(obs.hasNext()) { 
			sum += logLikelihood(obs.next());
		}
		return sum;
	}
	
	public double logLikelihood(Model m) { 
		BNValues parvals = new BNValues(m, parents);
		Integer childValue = child.encode(child.findValue(m));
		return Math.log(cpd.get(parvals).getProb(childValue));
	}
	
	public int countParameters() { 
		int count = child.size();
		for(int i = 0; i < parents.length; i++) { 
			count *= parents[i].size();
		}
		return count;
	}
	
	public Object sample(BNValues parvalues) { 
		return child.decode(cpd.get(parvalues).sampleValue());
	}
	
	public void resample(Model m) { 
		Object value = sample(new BNValues(m, parents));
		child.setValue(m, value);
	}
	
	public void learn(Iterator<? extends Model> obs) {
		for(BNValues vals : cpd.keySet()) { 
			cpd.get(vals).clear();
		}
		
		while(obs.hasNext()) { 
			Model m = obs.next();
			BNValues parvals = new BNValues(m, parents);
			Integer childValue = child.encode(child.findValue(m));
			cpd.get(parvals).addValue(childValue);
		}

		for(BNValues vals : cpd.keySet()) { 
			cpd.get(vals).normalize();
		}
	}
}
