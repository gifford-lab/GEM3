package edu.mit.csail.cgs.clustering.affinitypropagation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Vector;

public class CorrelationSimilarity {
	
	public static Map<String, List<Double>> fromFile(String infile) {
		Map<String, List<Double>> toreturn = new HashMap<String, List<Double>>();
		String line = "init";
		String[] splitline = new String[1];
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(infile));
			while (!((line = dataIn.readLine())==null)) {
				splitline = line.split("\t");
				List<Double> tmplist = new ArrayList<Double>();
				for (int i=1; i<splitline.length; i++) {
					tmplist.add(Double.valueOf(splitline[i]));
				}
				toreturn.put(splitline[0], tmplist);
			}
			dataIn.close();
		} catch (Exception e) {
			System.out.println("EXCEPTION");
		}
		System.out.println("Map size: "+toreturn.size());
		return toreturn;
	}
	
	public static void toFile(Map<String, List<Double>> valuemap, String outfile) {
		List<String> names = new ArrayList<String>(valuemap.keySet());
		double percent = 0;
		long[] counts = new long[5];
		for (int i=0; i<counts.length; i++) {
			counts[i] = 0;
		}
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			//PrintWriter out2 = new PrintWriter(new BufferedWriter(new FileWriter(outfile+".reduced")));
			System.out.println("Number of Names: "+names.size());
			for (int i=0; i<names.size(); i++) {
				for (int j=i+1; j<names.size(); j++) {
					if (i==j) continue;
					out.print(names.get(i));
					out.print("  ");
					out.print(names.get(j));
					out.print("  ");
					//out2.print(names.get(i));
					//out2.print("  ");
					//out2.print(names.get(j));
					//out2.print("  ");
					double value = convertCorrelationToSimilarity(computeSimilarity2(valuemap.get(names.get(i)), valuemap.get(names.get(j))));
					out.print(value);
					out.println();
					if (value > -0.1) {
						counts[0]++;
					}
					if (value > -0.2) {
						counts[1]++;
					}
					if (value > -0.3) {
						counts[2]++;
					}
					if (value > -0.4) {
						counts[3]++;
					}
					if (value > -0.5) {
						counts[4]++;
						//out2.print(value);
						//out2.println();
					}
				}
				if (i%(names.size() / 100)==0) {
					System.out.println(percent+"% done");
					percent += 1;
				}
			}
			out.flush();
			out.close();
			//out2.flush();
			//out2.close();
			for (int i=0; i<counts.length; i++) {
				System.out.println(counts[i]);
			}
			System.out.println();
		} catch (IOException e) {
			System.out.println("EXEPTION");
			e.printStackTrace();
		}
	}
	
	/**
	 * This method of computing correlation doesn't seem to work
	 * @param x
	 * @param y
	 * @return
	 */
	public static double computeSimilarity(List<Double> x, List<Double> y) {
		printlist(x);
		printlist(y);
		
		double sum_sq_x = 0.0;
		double sum_sq_y = 0.0;
		double sum_coproduct = 0.0;
		double mean_x = x.get(0);
		double mean_y = y.get(0);
		for (int i=1; i<x.size(); i++) {
			double sweep = (((double) i) - 1.0) / ((double) i);
			double delta_x = x.get(i) - mean_x;
			double delta_y = y.get(i) - mean_y;
			sum_sq_x += delta_x * delta_x * sweep;
			sum_sq_y += delta_y * delta_y * sweep;
			sum_coproduct += delta_x * delta_y * sweep;
			mean_x += delta_x / ((double)i);
			mean_y += delta_y / ((double)i);
		}
		double pop_sd_x = Math.sqrt( sum_sq_x / ((double)x.size()) );
		double pop_sd_y = Math.sqrt( sum_sq_y / ((double)x.size()) );
		double cov_x_y = sum_coproduct / ((double)x.size());
		double correlation = cov_x_y / (pop_sd_x * pop_sd_y);
		return correlation;
	}
	
	/**
	 * This seems to be the better way of computing correlation
	 * @param x
	 * @param y
	 * @return
	 */
	public static double computeSimilarity2(List<Double> x, List<Double> y) {
		int indx, n;
		double x_sum, y_sum, xx_sum, yy_sum, xy_sum;
		double[] dataPoint = new double[2];
		int X=0, Y=1;

		x_sum = y_sum = xx_sum = yy_sum = xy_sum = 0.0;
		indx = n = x.size();
		while (indx-- != 0) {
		    dataPoint[X] = x.get(indx);
		    dataPoint[Y] = y.get(indx);
		    x_sum  += dataPoint[X];
		    y_sum  += dataPoint[Y];
		    xx_sum += dataPoint[X] * dataPoint[X];
		    yy_sum += dataPoint[Y] * dataPoint[Y];
		    xy_sum += dataPoint[X] * dataPoint[Y];
		}

		return ((n * xy_sum - x_sum * y_sum) /
		     (Math.sqrt((n * xx_sum - x_sum * x_sum) *
				(n * yy_sum - y_sum * y_sum))));
	}
	
	public static void convertFile(String infile, String outfile) {
		String line = "init";
		String[] splitline = new String[1];
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(infile));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			while (!((line = dataIn.readLine())==null)) {
				splitline = line.split("  ");
				out.print(splitline[0]);
				out.print("  ");
				out.print(splitline[1]);
				out.print("  ");
				out.println(convertCorrelationToSimilarity(Double.valueOf(splitline[2])));
			}
			out.flush();
			out.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void reduceFile(String infile, String outfile) {
		long oldlines = 0;
		long newlines = 0;
		int perc = 0;
		String line = "init";
		String[] splitline = new String[1];
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(infile));
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			while (!((line = dataIn.readLine())==null)) {
				oldlines++;
				splitline = line.split("  ");
				if (Double.valueOf(splitline[2]).doubleValue() > -0.1) {
					out.println(line);
					newlines++;
				}
				if (oldlines%2E7==0) {
					System.out.println(perc+"% done");
					perc++;
				}
			}
			System.out.println(oldlines+" lines in old file");
			System.out.println(newlines+" lines in new file");
			out.flush();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static double convertCorrelationToSimilarity(double c) {
		return -(1.0 - c) / 2.0;
	}
	
	public static void generatePrefFile2(String infile, String outfile) {
		String line = "init";
		String[] splitline = new String[1];
		double[] values = new double[1];
		double median = 0.0;
		Vector<Double> valuevec = new Vector<Double>(81596385);
		int idx = 0;
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(infile));
			while (!((line = dataIn.readLine())==null)) {
				splitline = line.split("  ");
				valuevec.add(Double.valueOf(splitline[2]));
			}
			values = toDoubleArray(valuevec);
			Arrays.sort(values);
			if (values.length % 2 == 0) {
				median =  (values[values.length/2-1] + values[values.length/2]) / 2.0;
			} else {
				median =  values[(int)Math.floor(((double)values.length)/2.0)];
			}
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			for (int i=0; i<45101; i++) {
				out.println(median);
			}
			out.flush();
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void generatePrefFile3(double value, String outfile) {
		try {
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			for (int i=0; i<45101; i++) {
				out.println(value);
			}
			out.flush();
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static double[] toDoubleArray(Vector<Double> vec) {
		double[] toreturn = new double[vec.size()];
		for (int i=0; i<toreturn.length; i++) {
			toreturn[i] = vec.get(i);
		}
		return toreturn;
	}
	
	public static void generatePrefFile(String infile, String outfile) {
		String line = "init";
		String[] splitline = new String[1];
		TreeSet<Double> ts = new TreeSet<Double>();
		double median = 0.0;
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(infile));
			while (!((line = dataIn.readLine())==null)) {
				splitline = line.split("  ");
				ts.add(Double.valueOf(splitline[2]));
			}
			dataIn.close();
			Iterator<Double> it = ts.iterator();
			int half = (ts.size() / 2) - 1;
			for (int i=0; i < half; i++) {
				it.next();
			}
			if (ts.size()%2==0) {
				median = (it.next() + it.next()) / 2.0; 
			} else {
				it.next();
				median = it.next();
			}
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(outfile)));
			for (int i=0; i<45101; i++) {
				out.println(median);
			}
			out.flush();
			out.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void printlist(List l) {
		for (int i=0; i<l.size(); i++) {
			System.out.print(l.get(i));
			System.out.print(" ");
		}
		System.out.println();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		/*
		toFile(fromFile("/Users/reeder/tmp/lineage430ab_GCRMA_nohead.txt"), "/Users/reeder/tmp/Complete_Similarities.txt");
		*/
		/*
		convertFile("/afs/csail.mit.edu/u/c/ccr/psrg/projects/Corinna_Macklis/Data/lineage430ab_GCRMA_similarities.txt", "/Users/reeder/Corinna_similarities.txt");
		*/
		/*
		generatePrefFile2("/Users/reeder/Macklis_Proj/Complete_sims_reduced2.txt", "/Users/reeder/Macklis_Proj/Complete_prefs_red2.txt");
		*/
		/*
		reduceFile("/Users/reeder/Macklis_Proj/Complete_similarities.txt", "/Users/reeder/Macklis_Proj/Complete_sims_reduced2.txt");
		*/
		generatePrefFile3(-0.05, "/Users/reeder/Macklis_Proj/Complete_prefs_red2.txt");
	}

}
