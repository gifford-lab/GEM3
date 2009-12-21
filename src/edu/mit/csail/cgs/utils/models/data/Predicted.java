/*
 * Author: tdanford
 * Date: Aug 27, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import java.util.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.BitVector;
import edu.mit.csail.cgs.utils.models.*;
import Jama.*;

public class Predicted<M extends Model> {

	private DataFrame<M> frame;
	private Field numericField;
	
	public Predicted(DataFrame<M> f, String pfn) { 
		frame = f;
		numericField = null;
		
		Class<M> mcls = frame.getModelClass();

		try {
			Field mf = mcls.getField(pfn);
			Class t = mf.getType();
			if(Model.isSubclass(t, Number.class)) { 
				numericField = mf;
			} else { 
				throw new IllegalArgumentException(String.format(
						"Field %s is not numeric", pfn, t.getName()));					
			}
			
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(String.format(
					"Unknown field name: %s", pfn));					
		}
	}
	
	public Matrix createVector(BitVector selector) {
		int rows = selector != null ? selector.countOnBits() : frame.size();
		Matrix m = new Matrix(rows, 1);

		for(int i = 0, j = 0; j < frame.size(); j++) {
			if(selector == null || selector.isOn(j)) { 
				M rowValue = frame.object(j);
				
				try {
					Number numValue = (Number)numericField.get(rowValue);
					m.set(i, 0, numValue.doubleValue());

				} catch (IllegalAccessException e) {
					e.printStackTrace();
					throw new IllegalStateException(String.format("Couldn't access field %s: %s", 
							numericField.getName(), e.getMessage()));
				}
				
				i++;
			}
		}

		return m;
	}
	
	public Matrix createVector() { 
		int rows = frame.size();
		Matrix m = new Matrix(rows, 1);

		for(int i = 0; i < rows; i++) { 
			M rowValue = frame.object(i);
			try {
				Number numValue = (Number)numericField.get(rowValue);
				m.set(i, 0, numValue.doubleValue());

			} catch (IllegalAccessException e) {
				e.printStackTrace();
				throw new IllegalStateException(String.format("Couldn't access field %s: %s", 
						numericField.getName(), e.getMessage()));
			}
		}

		return m;
	}
}
