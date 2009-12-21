/*
 * Author: tdanford
 * Date: Nov 24, 2008
 */
package edu.mit.csail.cgs.utils.numeric;

import java.io.*;

import Jama.*;

public class JamaUtils {
	
	/*
	 * I/O Methods
	 */

	public static void printMatrix(Matrix X) { 
		printMatrix(X, System.out);
	}
	
	public static void printMatrix(Matrix X, PrintStream ps) { 
		printMatrix(X, ps, 2);
	}
	
	public static void printMatrix(Matrix X, PrintStream ps, int precision) { 
		for(int i = 0; i < X.getRowDimension(); i++) { 
			for(int j = 0; j < X.getColumnDimension(); j++) { 
				ps.print(String.format(String.format("%%.%df ", precision), X.get(i, j)));
			}
			ps.println();
		}
	}
	
	/**
	 * Normalizes a matrix by re-centering and re-scaling each of its columns, 
	 * so that they have a mean of 0 and a standard deviation of 1.
	 * 
	 * @param X
	 */
	public static void standardizeMatrix(Matrix X) { 
		Matrix mean = new Matrix(1, X.getColumnDimension(), 0.0);

		for(int i = 0; i < X.getRowDimension(); i++) { 
			for(int j = 0; j < X.getColumnDimension(); j++) { 
				mean.set(0, j, mean.get(0, j) + X.get(i, j));
			}
		}
		
		for(int j = 0; j < X.getColumnDimension(); j++) { 
			mean.set(0, j, mean.get(0, j) / (double)X.getRowDimension());
		}
		
		for(int i = 0; i < X.getRowDimension(); i++) { 
			for(int j = 0; j < X.getColumnDimension(); j++) { 
				X.set(i, j, X.get(i, j)-mean.get(0, j));
			}
		}
		
		Matrix std = new Matrix(1, X.getColumnDimension(), 0.0);

		for(int i = 0; i < X.getRowDimension(); i++) { 
			for(int j = 0; j < X.getColumnDimension(); j++) {
				double diff = X.get(i, j); 
				diff *= diff;
				std.set(0, j, std.get(0, j) + diff);
			}
		}
		
		for(int j = 0; j < X.getColumnDimension(); j++) { 
			std.set(0, j, Math.sqrt(std.get(0, j) / (double)X.getRowDimension()));
		}
		
		for(int i = 0; i < X.getRowDimension(); i++) { 
			for(int j = 0; j < X.getColumnDimension(); j++) { 
				X.set(i, j, X.get(i, j) / std.get(0, j));
			}
		}
	}
}
