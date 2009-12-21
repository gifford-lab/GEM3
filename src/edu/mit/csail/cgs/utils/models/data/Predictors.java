/*
 * Author: tdanford
 * Date: Aug 27, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import java.util.*;
import java.util.regex.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.BitVector;
import edu.mit.csail.cgs.utils.models.*;
import Jama.*;

public class Predictors<M extends Model> {

	private static Pattern interactionPattern = Pattern.compile("([^:]+):(.+)");
	
	private DataFrame<M> frame;
	private boolean hasConstant;
	
	private Vector<Field> numeric;  // the list of 'numeric' variables.
	
	private Vector<Field> factor;  // the list of factor variables.
	private Map<Field,Vector<Object>> factorCodes;  // the unique values that can be taken by each factor variable.
	
	// Each interaction is the product of one or more variables, some of which may be factors and the others 'numeric'.
	private Vector<Interaction> interactions;
	
	// Each InteractionValue is a product of values from
	// the variables in an interaction which are 'factor' variables.  
	private Map<Interaction,Vector<InteractionValue>> interactionValues;  
	
	private int cols;
	private Vector<String> columnNames;
	
	public Predictors(DataFrame<M> f, String... fs) { 
		frame = f;
		numeric = new Vector<Field>();
		factor = new Vector<Field>();
		interactions = new Vector<Interaction>();
		hasConstant = false;
		factorCodes = new HashMap<Field,Vector<Object>>();
		interactionValues = new HashMap<Interaction,Vector<InteractionValue>>();
		columnNames = new Vector<String>();
		cols = 0;
		
		Class<M> mcls = frame.getModelClass();
		Set<String> seenFields = new HashSet<String>();
		ModelFieldAnalysis<M> analysis = new ModelFieldAnalysis<M>(f.getModelClass());
		
		for(int i = 0; i < fs.length; i++) { 
			if(seenFields.contains(fs[i])) { 
				throw new IllegalArgumentException(String.format(
						"Duplicate field name: %s", fs[i]));
			}
			
			if(fs[i].equals("1")) { 
				hasConstant = true;
				cols += 1;
			} else { 
				try { 
					Matcher inMatcher = interactionPattern.matcher(fs[i]);
					Field mf = analysis.findField(fs[i]);
					if(mf != null) { 
						Class t = mf.getType();
						if(Model.isSubclass(t, Number.class)) { 
							addPredictor(fs[i]);
						} else if (Model.isSubclass(t, String.class)) { 
							addFactor(fs[i]);
						} else { 
							throw new IllegalArgumentException(String.format(
									"Field %s is not a regression-ready predictor", fs[i]));
						}

					} else if (inMatcher.matches()) { 
						String[] array = fs[i].split(":");
						addInteraction(array);
					} else { 
						throw new IllegalArgumentException(String.format(
								"Unknown field name: %s", fs[i]));						
					}
					
				} catch(NoSuchFieldException e) { 
					throw new IllegalArgumentException(String.format(
							"Unknown field name: %s", fs[i]));
				}
			}
			
			seenFields.add(fs[i]);
		}
		
		if(hasConstant) { 
			columnNames.insertElementAt("(Intercept)", 0);
		}
	}
	
	public int size() {
		return frame.size();
	}

	public void addConstant() {
		if(!hasConstant) { 
			hasConstant = true;
			columnNames.insertElementAt("(Intercept)", 0);
		}
	}
	
	public void addInteraction(String... fns) throws NoSuchFieldException {
		ModelFieldAnalysis<M> mfa = new ModelFieldAnalysis<M>(frame.getModelClass());
		Field[] fields = mfa.findFields(fns);
		for(int i = 0; i < fields.length; i++) { 
			if(fields[i] == null) { throw new NoSuchFieldException(fns[i]); }
		}
		
		Interaction inter = new Interaction(fields);
		
		if(interactions.contains(inter)) { 
			throw new IllegalArgumentException(
					String.format("Cannot add the same interaction %s twice.", inter.toString()));
		}
		interactions.add(inter);
		
		String[] factorFields = inter.findFactorFieldNames();
		Vector<InteractionValue> values = inter.allInteractionValues(factorFields);
		
		interactionValues.put(inter, values);
		cols += values.size();
		
		for(InteractionValue value : values) { 
			String colName = String.format("%s(%s)", inter.toString(), value.toString());
			columnNames.add(colName);
		}
	}
	
	public void addPredictor(String fn) throws NoSuchFieldException { 
		Field f = frame.getModelClass().getField(fn);
		Class t = f.getType();
		if(Model.isSubclass(t, Double.class)) { 
			numeric.add(f);
		} else if (Model.isSubclass(t, Integer.class)) { 
			numeric.add(f);
		} else { 
			throw new NoSuchFieldException(String.format(
					"%s is not a numeric field (%s)", fn, t.getName()));
		}
		cols += 1;
		columnNames.add(fn);
	}
	
	public Set<String> findFactorValues(String fn) { 
		TreeSet<String> values = new TreeSet<String>();
		for(Object v : frame.fieldValues(fn)) { 
			values.add((String)v);
		}
		values.remove(values.first());
		return values;
	}
	
	public void addFactor(String fn) throws NoSuchFieldException {
		Field f = frame.getModelClass().getField(fn);
		Class t = f.getType();
		if(Model.isSubclass(t, String.class)) {
			Set<String> values = findFactorValues(fn);
			factor.add(f);			
			//factorCodes.put(f, new Vector(new TreeSet<String>(values)));
			factorCodes.put(f, new Vector(values));
			cols += values.size();
			for(String obj : values) { 
				columnNames.add(String.format("%s(%s)", fn, obj));
			}
			
		} else { 
			throw new IllegalArgumentException(String.format("%s is not a valid factor-field.", fn));
		}
	}
	
	public int getNumColumns() { return cols; }
	public String getColumnName(int i) { return columnNames.get(i); }
	
	public Matrix createMatrix() { 
		return createMatrix(null);
	}
	
	public Matrix createMatrix(BitVector selector) { 
		return createMatrix(selector, null);
	}
	
	public Matrix createMatrix(BitVector selector, Map<String,Transformation<Double,Double>> transforms) { 
		int rows = selector != null ? selector.countOnBits() : frame.size();
		Matrix m = new Matrix(rows, cols);
		
		int cidx = 0;
		if(hasConstant) { 
			for(int i = 0, j = 0; j < frame.size(); j++) {
				if(selector == null || selector.isOn(j)) { 
					m.set(i, cidx, 1.0);
					i += 1;
				}
			}
			cidx += 1;
		}

		for(Field f : numeric) {
			Transformation<Double,Double> transform = 
				transforms != null && transforms.containsKey(f.getName()) ? 
					transforms.get(f.getName()) : null;
			
			for(int i = 0, j = 0; j < frame.size(); j++) {
				if(selector == null || selector.isOn(j)) { 
					M rowValue = frame.object(j);
					try {
						Double numValue = ((Number)f.get(rowValue)).doubleValue();
						if(transform != null) { 
							numValue = transform.transform(numValue);
						}
						
						m.set(i, cidx, numValue);

					} catch (IllegalAccessException e) {
						e.printStackTrace();
						throw new IllegalStateException(String.format("Couldn't access field %s: %s", 
								f.getName(), e.getMessage()));
					}
					i+=1;
				}
			}
			cidx += 1;
		}
		
		for(Field f : factor) {
			Vector<Object> factorValues = factorCodes.get(f);
			for(int i = 0, j = 0; j < frame.size(); j++) {
				if(selector==null || selector.isOn(j)) { 
					M rowValue = frame.object(j);
					try {
						Object factorValue = f.get(rowValue);
						int idx = factorValues.indexOf(factorValue);
						if(idx != -1) { 
							m.set(i, cidx+idx, 1.0);
						}

					} catch (IllegalAccessException e) {
						e.printStackTrace();
						throw new IllegalStateException(String.format("Couldn't access field %s: %s", 
								f.getName(), e.getMessage()));
					}
					i += 1;
				}
			}
			
			cidx += factorValues.size();
		}

		for(Interaction in : interactions) { 
			Vector<InteractionValue> values = interactionValues.get(in);
			for(int i = 0, j = 0; j < frame.size(); j++) {
				if(selector==null || selector.isOn(j)) { 
					M rowValue = frame.object(j);
					InteractionValue inValue = in.calculateInteractionValue(rowValue);
					Double numValue = in.calculatePredictor(rowValue);

					int idx = values.indexOf(inValue);
					if(idx != -1) { 
						m.set(i, cidx+idx, numValue);
					}
					i+=1;
				}
			}
			cidx += values.size();
		}
		
		return m;
	}
	
	private class Interaction { 
		
		public Set<Field> fields;
		
		public Interaction(Field... fs) { 
			fields = new HashSet<Field>();
			for(int i = 0; i < fs.length; i++) { 
				fields.add(fs[i]);
			}
		}
		
		public int hashCode() { 
			int code = 17;
			for(Field f : fields) { 
				code += f.hashCode();
			}
			code *= 37;
			return code;
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder();

			for(Field f : fields) { 
				if(isFactorField(f)) { 
					if(sb.length() > 0) { sb.append(":"); }
					sb.append(f.getName());
				}
			}
			
			for(Field f : fields) { 
				if(!isFactorField(f)) { 
					if(sb.length() > 0) { sb.append(":"); }
					sb.append(f.getName());
				}
			}
			
			return sb.toString();
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof Predictors.Interaction)) { return false; }
			Interaction in = (Interaction)o;
			if(fields.size() != in.fields.size()) { return false; }
			for(Field f : fields) { 
				if(!(in.fields.contains(f))) { 
					return false; 
				}
			}
			return true;
		}
		
		private boolean isNumericField(Field f) { 
			return Model.isSubclass(f.getType(), Number.class);
		}
		
		private boolean isFactorField(Field f) { 
			return Model.isSubclass(f.getType(), String.class);
		}
		
		public InteractionValue calculateInteractionValue(Object o) { 
			InteractionValue value = new InteractionValue();
			
			for(Field f : fields) { 
				if(isFactorField(f)) {  
					try {
						String v = (String) f.get(o);
						value.values.add(v);
						
					} catch (IllegalAccessException e) {
						e.printStackTrace();
					}
				}
			}
			return value;			
		}
		
		public Double calculatePredictor(Object o) {
			Double value = 1.0;
			for(Field f : fields) { 
				if(isNumericField(f)) {  
					try {
						Number n = (Number)f.get(o);
						value *= n.doubleValue();
					} catch (IllegalAccessException e) {
						e.printStackTrace();
					}
				}
			}
			return value;
		}
		
		public String[] findFactorFieldNames() { 
			Vector<String> fs = new Vector<String>();
			for(Field f : fields) { 
				if(isFactorField(f)) {  
					fs.add(f.getName());
				}
			}
			return fs.toArray(new String[fs.size()]);
		}
		
		public String[] findNumericFieldNames() { 
			Vector<String> fs = new Vector<String>();
			for(Field f : fields) { 
				if(isNumericField(f)) {  
					fs.add(f.getName());
				}
			}
			return fs.toArray(new String[fs.size()]);
		}
		
		public Vector<InteractionValue> allInteractionValues(String[] names) { 
			Vector<InteractionValue> vv = new Vector<InteractionValue>();
			vv.add(new InteractionValue());
			
			String[] ff = findFactorFieldNames();
			for(int i = 0; i < ff.length; i++) {
				//Set values = findFactorValues(ff[i]);
				Set values = frame.fieldValues(ff[i]);
				
				vv = appendValues(vv, values);
			}
			
			// Analogous to taking out the very first value of a set of factor values.
			vv.remove(0);
			
			return vv;
		}
		
		private Vector<InteractionValue> appendValues(Vector<InteractionValue> prev, Set vals) { 
			Vector<InteractionValue> newv = new Vector<InteractionValue>();
			for(InteractionValue v : prev) { 
				newv.addAll(v.extend(vals));
			}
			return newv;
		}
	}
	
	private class InteractionValue { 
		
		public Vector values;
		
		public InteractionValue() { 
			values = new Vector();
		}
		
		public InteractionValue(InteractionValue v ) { 
			values = new Vector(v.values);
		}
		
		public InteractionValue(InteractionValue v, Object o) { 
			values = new Vector(v.values);
			values.add(o);
		}
		
		public Vector<InteractionValue> extend(Set vals) { 
			Vector<InteractionValue> ivs = new Vector<InteractionValue>();
			for(Object v : vals) { 
				ivs.add(new InteractionValue(this, v));
			}
			return ivs;
		}
		
		public String toString() { 
			StringBuilder sb = new StringBuilder();
			for(int i = 0; i < values.size(); i++) { 
				if(i > 0) { sb.append("_"); }
				sb.append(values.get(i).toString());
			}
			return sb.toString();
		}
		
		public int hashCode() { 
			 int code = 17;
			 for(Object v : values) { code += v.hashCode(); code *= 37; }
			 return code;
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof Predictors.InteractionValue)) { return false; }
			InteractionValue iv = (InteractionValue)o;
			if(iv.values.size() != values.size()) { return false; }
			for(int i = 0; i < values.size(); i++) { 
				if(!values.get(i).equals(iv.values.get(i))) { return false; }
			}
			return true;
		}
	}

	public Vector<String> getColumnNames() {
		return columnNames;
	}
}


