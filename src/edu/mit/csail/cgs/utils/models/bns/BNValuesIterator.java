/*
 * Author: tdanford
 * Date: Dec 3, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import java.util.Iterator;

public class BNValuesIterator implements Iterator<BNValues> {
	
	private BNVar[] vars;
	private Integer[] encodedValues;
	private Object[] values;
	private BNValues nextValues;
	
	public BNValuesIterator(BNVar... vars) { 
		this.vars = vars.clone();
		encodedValues = null;
		values = new Object[vars.length];
		findNextValues();
	}
	
	public boolean hasNext() {
		return nextValues != null;
	}
	
	private void findNextValues() { 
		nextValues = null;
		
		if(encodedValues == null) { 
			encodedValues = new Integer[vars.length];
			for(int i = 0; i < encodedValues.length; i++) { 
				encodedValues[i] = 0;
			}
		} else { 
			int i = encodedValues.length-1;
			while(i >= 0 && encodedValues[i] == vars[i].size()-1) { 
				i--;
			}
			if(i >= 0) { 
				encodedValues[i] += 1;
				for(int j = i + 1; j < encodedValues.length; j++) { 
					encodedValues[j] = 0;
				}
			} else { 
				return;
			}
		}

		for(int i = 0; i < encodedValues.length; i++) { 
			values[i] = vars[i].decode(encodedValues[i]);
		}
		
		nextValues = new BNValues(vars, values);
	}

	public BNValues next() {
		BNValues vals = nextValues;
		findNextValues();
		return vals;
	}

	public void remove() {
		throw new UnsupportedOperationException("remove()");
	}

}
