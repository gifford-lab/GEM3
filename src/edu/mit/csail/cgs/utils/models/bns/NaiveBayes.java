/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.mit.csail.cgs.utils.models.bns;

import edu.mit.csail.cgs.utils.ArrayUtils;
import edu.mit.csail.cgs.utils.models.*;
import edu.mit.csail.cgs.utils.models.data.DataFrame;

public class NaiveBayes<X extends Model> extends BN<X> {

	public NaiveBayes(DataFrame<X> data, String classField, String... attrs) { 
		super(data, ArrayUtils.prepend(classField, attrs));
		for(int i = 0; i < attrs.length; i++) { 
			graph.addEdge(classField, attrs[i]);
		}
		learnCPDs();
	}
	
	
}
