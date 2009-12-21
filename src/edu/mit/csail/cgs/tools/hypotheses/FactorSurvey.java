package edu.mit.csail.cgs.tools.hypotheses;

import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.SetTools;

public class FactorSurvey {

	private Factor factor;
	private BindingExplorer explorer;
	
	public FactorSurvey(Factor f, BindingExplorer e) { 
		factor = f;
		explorer = e;
	}
	
	public Factor getFactor() { 
		return factor;
	}
	
	public HypothesisTree createTree(int depth) { 
		LinkedList<BindingHypothesis> hyps = createHypotheses(true);
		HypothesisTree ht = new HypothesisTree(explorer);
        ht.addChildrenAtDepth(hyps, depth);
		return ht;
	}
	
	public BindingExplorer getExplorer() { return explorer; }
	
	public LinkedList<BindingHypothesis> createHypotheses(boolean forward) { 
		LinkedList<BindingHypothesis> hyps = new LinkedList<BindingHypothesis>();
		for(int i = 0; i < explorer.getNumFactors(); i++) { 
			Factor f = explorer.getFactor(i);
			if(!f.equals(factor)) { 
				BindingHypothesis hyp = forward ? 
					new BindingHypothesis.Conditional(factor, f) : 
					new BindingHypothesis.Conditional(f, factor);
				hyps.add(hyp);
			}
		}
		return hyps;
	}
	
}
