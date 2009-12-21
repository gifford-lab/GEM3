package edu.mit.csail.cgs.tools.hypotheses;

import java.util.*;
import java.util.regex.*;
import edu.mit.csail.cgs.datasets.general.*;

public interface BindingHypothesis {
	
	public boolean isSupportedBy(LocalBindingSummary lbs);
	public Set<Factor> getInvolvedFactors();
    public int getOptimizationClass();
	
	public static class SimpleBinding implements BindingHypothesis {
		
		private Factor factor;
		
		public SimpleBinding(Factor f) { factor = f; }

		public Set<Factor> getInvolvedFactors() {
			HashSet<Factor> fs = new HashSet<Factor>();
			fs.add(factor);
			return fs;
		}

		public boolean isSupportedBy(LocalBindingSummary lbs) {
			return lbs.getBound().contains(factor);
		} 
        
        public Factor getFactor() { return factor; }
        
        public int getOptimizationClass() { return 0; }
		
		public String toString() { return "[" + factor.getName() + "]"; }
		
		public int hashCode() { 
			int code = 17;
			code += factor.hashCode(); code *= 37;
			return code;
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof SimpleBinding)) { return false; }
			SimpleBinding sb = (SimpleBinding)o;
			return factor.equals(sb.factor);
		}
	}
	
	public static class Conditional implements BindingHypothesis {
		
		private BindingHypothesis antecedent, consequent;
        private int optClass;
		
		public Conditional(BindingHypothesis ifh, BindingHypothesis thenh) { 
			antecedent = ifh;
			consequent = thenh;
            optClass = 0;
		}
		
		public Conditional(Factor f1, Factor f2) {
			antecedent = new SimpleBinding(f1);
			consequent = new SimpleBinding(f2);
            optClass = 1;
		}
        
        public BindingHypothesis getAntecedent() { return antecedent; }
        public BindingHypothesis getConsequent() { return consequent; }
        
        public int getOptimizationClass() { return optClass; }

		public Set<Factor> getInvolvedFactors() {
			HashSet<Factor> fs = new HashSet<Factor>();
			fs.addAll(antecedent.getInvolvedFactors());
			fs.addAll(consequent.getInvolvedFactors());
			return fs;
		}

		public boolean isSupportedBy(LocalBindingSummary lbs) {
			return consequent.isSupportedBy(lbs) || !antecedent.isSupportedBy(lbs);
		} 
		
		public String toString() { 
			return "(" + antecedent + " -> " + consequent + ")"; 
		}
		
		public int hashCode() { 
			int code = 17;
			code += antecedent.hashCode(); code *= 37;
			code += consequent.hashCode(); code *= 37;
			return code;
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof Conditional)) { return false; }
			Conditional c = (Conditional)o;
			if(!antecedent.equals(c.antecedent)) { return false; }
			if(!consequent.equals(c.consequent)) { return false; }
			return true;
		}
	}
}
