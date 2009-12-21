/*
 * Author: tdanford
 * Date: Aug 28, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.models.IllegalModelException;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.ModelFieldAnalysis;


public class RegressionModel extends Model {
	
	public Field getDependentVariable() { 
		ModelFieldAnalysis mfa = new ModelFieldAnalysis(getClass());
		Vector<Field> fs = mfa.findTypedFields(DependentVariable.class);
		if(fs.size() == 1) { 
			return fs.get(0); 
		} else if (fs.size() == 0) { 
			throw new IllegalModelException(String.format(
					"Model %s has no dependent variable", getClass().getName()));
		} else { 
			StringBuilder sb = new StringBuilder();
			for(Field f : fs) {
				if(sb.length() > 0) { sb.append(" "); }
				sb.append(f.getName());
			}
			throw new IllegalModelException(String.format(
					"Model %s has more than one dependent variable: %s", 
					getClass().getName(), sb.toString()));
		}
	}
	
	
	public boolean hasInterceptVariable() { 
		ModelFieldAnalysis mfa = new ModelFieldAnalysis(getClass());
		Vector<Field> fs = mfa.findTypedFields(Intercept.class);
		if(fs.size() > 1) { 
			StringBuilder sb = new StringBuilder();
			for(Field f : fs) {
				if(sb.length() > 0) { sb.append(" "); }
				sb.append(f.getName());
			}
			throw new IllegalStateException(
					String.format("Model has more than one intercept term (%s)", sb.toString()));
		}
		return fs.size() == 1;
	}
	
	
	public Vector<Field> getIndependentVariables() { 
		ModelFieldAnalysis mfa = new ModelFieldAnalysis(getClass());
		Vector<Field> fs = new Vector<Field>();
		fs.addAll(mfa.findTypedFields(NumericVariable.class));
		fs.addAll(mfa.findTypedFields(FactorVariable.class));
		return fs;
	}
	
	public static interface DependentVariable {
	}
	
	public static interface FactorVariable { 
	}
	
	public static interface NumericVariable { 
	}
	
	public static interface Intercept { 
	}
}
