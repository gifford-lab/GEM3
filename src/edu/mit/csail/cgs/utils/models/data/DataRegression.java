/*
 * Author: tdanford
 * Date: Aug 27, 2008
 */
package edu.mit.csail.cgs.utils.models.data;

import java.util.*;
import java.io.*;
import java.util.regex.*;
import java.lang.reflect.*;

import cern.jet.random.ChiSquare;
import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

import edu.mit.csail.cgs.utils.BitVector;
import edu.mit.csail.cgs.utils.Predicate;
import edu.mit.csail.cgs.utils.models.*;
import Jama.*;

public class DataRegression<M extends Model> {
	
	public static void main(String[] args) { 
		File f = new File("C:\\Documents and Settings\\tdanford\\Desktop\\test.txt");
		try {
			DataFrame<XYPoint> df = new DataFrame<XYPoint>(XYPoint.class, f);

			RegressionModel m = new RegressionModel() { 
				public DependentVariable y;
				public NumericVariable x;
				public Intercept b;
			};
			
			DataRegression<XYPoint> reg = new DataRegression<XYPoint>(df, "y ~ x + 1");
			//DataRegression<XYPoint> reg = new DataRegression<XYPoint>(df, m);
			
			reg.transform(new ATransformation<XYPoint,XYPoint>(XYPoint.class,XYPoint.class) {
				public XYPoint transform(XYPoint v) {
					//v.x -= 1.0;
					v.y *= 2.0;
					return v;
				} 
			});
			
			Map<String,Double> coeffs = reg.calculateRegression();
			Map<String,Double[]> bounds = reg.calculateBounds();
			
			for(String title : coeffs.keySet()) {
				Double[] b = bounds.get(title);
				System.out.println(String.format("%s \t%.3f\t(%.3f, %.3f)", 
						title, coeffs.get(title), 
						b[0], b[1]));
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private DataFrame<M> frame;
	private Predicted<M> dataY;
	private Predictors<M> dataX;
	
	private QRDecomposition qr;
	private Matrix Rinv, betaHat, Vbeta;
	private double s2, r2;
	
	private String dataYVar;
	private String[] dataXVars;
	
	private RandomEngine engine;
	private Normal ndist;
	
	public DataRegression(DataFrame<M> f, String stmt) { 
		frame = f;
		engine = new cern.jet.random.engine.DRand();
		ndist = new Normal(0.0, 1.0, engine);
		
		Vector<String> vs = parseStatement(stmt);
		if(vs == null) { 
			throw new IllegalArgumentException(String.format("Couldn't parse statement \"%s\"",
					stmt));
		}
		
		dataYVar = vs.get(0);
		dataXVars = vs.subList(1, vs.size()).toArray(new String[vs.size()-1]);
		
		dataY = new Predicted<M>(frame, dataYVar);
		dataX = new Predictors<M>(frame, dataXVars);
	}
	
	public DataRegression(DataFrame<M> f, RegressionModel m) { 
		frame = f;
		engine = new cern.jet.random.engine.DRand(); 
		ndist = new Normal(0.0, 1.0, engine);
		
		Field dvar = m.getDependentVariable();
		Vector<Field> ivars = m.getIndependentVariables();
		boolean intercept = m.hasInterceptVariable();
		int plus = intercept ? 1 : 0;
		
		dataYVar = dvar.getName();
		dataXVars = new String[ivars.size() + plus];
		
		int i = 0;
		if(intercept) { dataXVars[i++] = "1"; }
		for(; i < dataXVars.length; i++) { 
			dataXVars[i] = ivars.get(i-plus).getName();
		}
		
		dataY = new Predicted<M>(frame, dataYVar);
		dataX = new Predictors<M>(frame, dataXVars);
	}
	
	public void filter(Predicate<M> p) { 
		frame.filter(p);
	}
	
	public void transform(Transformation<M,M> t) {
		frame = frame.transform(t);
		dataY = new Predicted<M>(frame, dataYVar);
		dataX = new Predictors<M>(frame, dataXVars);		
	}
	
	public Vector<String> getPredictorNames() { 
		return dataX.getColumnNames();
	}
	
	public Predictors<M> getPredictors() { return dataX; }
	public Predicted<M> getPredicted() { return dataY; }
	
	public Matrix getPredictorMatrix() { 
		return dataX.createMatrix();
	}
	
	public Matrix getPredictedVector() { 
		return dataY.createVector();
	}
	
	public DataFrame<M> getFrame() { return frame; }
	
	public Map<String,Double> calculateRegression() { 
		calculate();
		return collectCoefficients();
	}
	
	public Map<String,Double> collectCoefficients() { 
		HashMap<String,Double> map = new LinkedHashMap<String,Double>();
		for(int i = 0; i < betaHat.getRowDimension(); i++) { 
			String name = dataX.getColumnName(i);
			map.put(name, betaHat.get(i, 0));
		}
		
		return map;		
	}
	
	public Map<String,Double[]> calculateBounds() { 
		HashMap<String,Double[]> map = new LinkedHashMap<String,Double[]>();
		Vector<Double[]> bounds = sampleBetaBounds(100);
		
		for(int i = 0; i < betaHat.getRowDimension(); i++) { 
			String name = dataX.getColumnName(i);
			map.put(name, bounds.get(i));
		}
		
		return map;
		
	}
	
	public void calculate() { 
		calculate(null);
	}
	
	public void calculate(BitVector selector) { 
		calculate(selector, null);
	}
	
	public void calculate(BitVector selector, Map<String,Transformation<Double,Double>> transforms) {
		Matrix y = dataY.createVector(selector);
		Matrix X = dataX.createMatrix(selector, transforms);
		calculate(X, y); 
	}
	
	public static Matrix leastSquares(Matrix X, Matrix y) { 
		QRDecomposition qr = new QRDecomposition(X);
		Matrix R = qr.getR();
		Matrix Qtransy = qr.getQ().transpose().times(y);
		Matrix betaHat = R.solve(Qtransy);
		
		return betaHat;
	}
	
	public static double s2(Matrix X, Matrix y, Matrix betaHat) { 
		Matrix yhat = X.times(betaHat);
		Matrix errors = y.minus(yhat);
		int n = X.getRowDimension(), k = X.getColumnDimension();
		
		double s2 = (errors.transpose().times(errors)).get(0, 0);
		s2 /= (double)(n - k);

		return s2;
	}
		
	public void calculate(Matrix X, Matrix y) { 
		qr = new QRDecomposition(X);
		Matrix R = qr.getR();
		Rinv = R.inverse();
		Vbeta = Rinv.times(Rinv.transpose());
		
		Matrix Qtransy = qr.getQ().transpose().times(y);
		betaHat = R.solve(Qtransy);
		
		Matrix yhat = X.times(betaHat);
		Matrix errors = y.minus(yhat);
		int n = X.getRowDimension(), k = X.getColumnDimension();
		
		s2 = (errors.transpose().times(errors)).get(0, 0);
		s2 /= (double)(n - k);
		
		calculateR2(y, yhat);
	}
	
	private void calculateR2(Matrix y, Matrix yhat) { 
		double mean = 0.0;
		
		for(int i = 0; i < y.getRowDimension(); i++) { 
			double yvalue = y.get(i, 0);
			mean += yvalue;
		}
		
		mean /= (double)y.getRowDimension();
		
		double SSE = 0.0, SST = 0.0, SSR = 0.0;
		for(int i = 0; i < y.getRowDimension(); i++) { 
			double yvalue = y.get(i, 0);
			double sstDiff = yvalue-mean;
			double sseDiff = yvalue-yhat.get(i, 0);
			double ssrDiff = yhat.get(i,0) - mean;
			
			SST += (sstDiff * sstDiff);
			SSE += (sseDiff * sseDiff);
			SSR += (ssrDiff * ssrDiff);
		}
		
		r2 = 1.0 - (SSE / SST);
	}
	
	public Matrix getBetaHat() { return betaHat; }
	public Matrix getVarBeta() { return Vbeta; }
	public double getR2() { return r2; } 
	public double getS2() { return s2; }
	public int getN() { return dataX.size(); }
	public int getK() { return dataX.getNumColumns(); }
	
	public Vector<Double[]> sampleBetaBounds(int iters) { 
		Vector<Double[]> v = new Vector<Double[]>();
		for(int i = 0; i < getK(); i++) { 
			v.add(new Double[iters]); 
		}
		
		for(int i = 0; i < iters; i++) {
			double var = sampleVar();
			Matrix beta = sampleBeta(var);
			
			for(int j = 0; j < getK(); j++) { 
				v.get(j)[i] = beta.get(j, 0);
			}
		}
		
		Vector<Double[]> bounds = new Vector<Double[]>();
		int lower = (iters/4);
		int upper = 3*(iters/4);
		
		for(int j = 0; j < getK(); j++) { 
			Double[] sarray = v.get(j);
			Arrays.sort(sarray);
			Double[] b = new Double[] { sarray[lower], sarray[upper] };
			bounds.add(b);
		}
		
		return bounds;
	}
	
	public Matrix sampleBeta(double var) {
		Matrix beta = new Matrix(getK(), 1);
		for(int i = 0; i < beta.getRowDimension(); i++) {
			double n = ndist.nextDouble();
			beta.set(i, 0, n);
		}
		double sd = Math.sqrt(var);
		beta = Rinv.times(sd).times(beta).plus(betaHat);
		return beta;
	}
	
	public double sampleVar() { 
		double diff = (double)(getN() - getK());
		ChiSquare chiSquare = new cern.jet.random.ChiSquare(diff, engine);
		double x = chiSquare.nextDouble();
		return (diff * s2) / x;
	}
	
	/**
	 * @deprecated
	 * @param betaHat
	 * @return
	 */
	public double calculateS2(Matrix betaHat) { 
		// pg. 356 of Gelman

		// n : number of datapoints
		// k : number of predictors
		double n = (double)frame.size();
		double k = (double)dataX.getNumColumns();
		double coeff = 1.0 / (n - k);
		
		Matrix y = dataY.createVector();
		Matrix X = dataX.createMatrix();
		if(betaHat == null) { betaHat = calculateBetaHat(); }
		
		Matrix half = y.minus(X.times(betaHat));
		
		Matrix product = half.transpose().times(half);
		
		double ret = coeff * product.get(0, 0);
		return ret;
	}
	
	/**
	 * @deprecated
	 * @return
	 */
	public Matrix calculateBetaHat() {
		// pg. 356 of Gelman
		
		// n : number of datapoints
		// k : number of predictors
		
		// Xtrans : k x n
		Matrix Xtrans = dataX.createMatrix().transpose();
		
		// Vbeta : k x k 
		Matrix Vbeta = Xtrans.times(Xtrans.transpose());
		Vbeta = Vbeta.inverse();
		
		// ytransf : k x 1 
		Matrix ytransf = Xtrans.times(dataY.createVector());
		
		// res : k x 1
		Matrix res = Vbeta.times(ytransf);
		
		return res;
	}

	private static Pattern stmtPattern = Pattern.compile(
	"\\s*([^\\s~]+)\\s*~\\s*(.*)");

	private static Vector<String> parseStatement(String stmt) {  
		Matcher m = stmtPattern.matcher(stmt);
		Vector<String> v = null; 
		if(m.matches()) { 
			v = new Vector<String>();
			String y = m.group(1);
			v.add(y);

			String preds = m.group(2);
			String[] array = preds.split("\\s+");
			boolean seenConstant = false;
			boolean omitConstant = false;
			
			boolean lastMinus = false;

			for(int i = 0; i < array.length; i++) {
				if(i % 2 == 1) { 
					if(array[i].equals("-")) { 
						lastMinus = true;
					} else if(array[i].equals("+")) { 
						lastMinus = false;
					} else { 
						return null;
					}
				} else {
					if(array[i].equals("1")) { 
						seenConstant = true;
						if(lastMinus) { 
							omitConstant = true;
						} else { 
							v.add(array[i]);
						}
					} else { 
						v.add(array[i]);
					}
				}
			}

			if(!seenConstant && !omitConstant) { 
				v.add("1");
			}
		}
		return v;
	}
	
	public static void printMatrix(Matrix m, PrintStream ps, int precision) {
		String format = "%." + precision + "f";
		
		ps.print("   \t");
		for(int j = 0; j < m.getColumnDimension(); j++) { 
			if(j > 0) { ps.print("  "); }
			ps.print(String.format(" %3d", j));
		}
		ps.println();
		
		for(int i = 0; i < m.getRowDimension(); i++) {
			ps.print(String.format("%3d\t", i));
			for(int j = 0; j < m.getColumnDimension(); j++) { 
				if(j > 0) { ps.print("  "); }
				ps.print(String.format(format, m.get(i, j)));
			}
			ps.println();
		}
	}

}

