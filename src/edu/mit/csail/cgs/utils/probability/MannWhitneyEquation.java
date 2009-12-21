/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.utils.probability;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 * 
 * This implementation is coded from the paper:
 * 
 * Mann, H.B., and Whitney, D.R.. "On a Test of Whether One 
 * of Two Random Variables is Stochastically Larger Than the Other."
 * The Annals of Mathematical Statistics. Vol. 18, No. 1.  
 * (March 1947), pp.50-60.
 * 
 * URL: http://links.jstor.org/sici?sici=0003-4851%28194703%2918%3A1%3C50%3AOATOWO%3E2.0.CO%3B2-M
 * 
 */
public class MannWhitneyEquation implements RecurrenceEquation {

	public static void main(String[] args) { 
		explicit_main(args);
		//bits_main(args);
	}

	public static void explicit_main(String[] a) { 
		String line;
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.print("> "); System.out.flush();
		MannWhitneyEquation eq = new MannWhitneyEquation();
		int[] argArray = new int[3];
		RecurrenceEquation.Arguments args = null;

		try {
			while((line = br.readLine()) != null) { 
				line = line.trim();
				if(line.length() > 0) { 
					String[] array = line.split("\\s+");
					int n = Integer.parseInt(array[0]);
					int m = Integer.parseInt(array[1]);
					int u = Integer.parseInt(array[2]);

					double pvalue = eq.getLowerPValue(m, n, u);

					System.out.println("N: " + n + ", M: " + m + ", U: " + u);
					System.out.println("Lower P-value: " + pvalue);
				}

				System.out.print("> "); System.out.flush();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void bits_main(String[] a) { 
		String line;
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.print("> "); System.out.flush();
		MannWhitneyEquation eq = new MannWhitneyEquation();
		int[] argArray = new int[3];
		RecurrenceEquation.Arguments args = null;

		try {
			while((line = br.readLine()) != null) { 
				line = line.trim();
				if(line.length() > 0) { 
					int n = calculateN(line);
					int m = calculateM(line);
					int u = calculateU(line);
					//int u = calculateUFromT(m, n, line);

					System.out.println("N: " + n + ", M: " + m + ", U: " + u);
					System.out.println("Prob: " + eq.getProb(m, n, u));
					System.out.println("Lower P-value: " + eq.getLowerPValue(m, n, u));
					System.out.println("Upper P-value: " + eq.getUpperPValue(m, n, u));
				}
				System.out.print("> "); System.out.flush();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private Map<Arguments,Double> cache;

	public MannWhitneyEquation() {
		super();
		cache = new HashMap<Arguments,Double>();
	}

	public MannWhitneyEquation(File f) {
		super();
		cache = parseCacheFromFile(f);
	}

	private Map<Arguments,Double> parseCacheFromFile(File f) {
		Map<Arguments,Double> tmpCache = new HashMap<Arguments,Double>();
		try {

			BufferedReader br = new BufferedReader(new FileReader(f));
			String line = null;
			String[] split;
			while((line = br.readLine()) != null) {
				split = line.split("\t");
				int[] tmpArray = new int[3];
				tmpArray[0] = Integer.parseInt(split[0]);
				tmpArray[1] = Integer.parseInt(split[1]);
				tmpArray[2] = Integer.parseInt(split[2]);
				Arguments tmpArgs = new Arguments(tmpArray);
				tmpCache.put(tmpArgs, Double.parseDouble(split[3]));
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
		return tmpCache;
	}

	public static int calculateM(String bits) { 
		int count = 0; 
		for(int i = 0; i < bits.length(); i++) { 
			if(bits.charAt(i) == '0') { 
				count += 1; 
			}
		}
		return count;
	}

	public static int calculateN(String bits) { 
		int count = 0; 
		for(int i = 0; i < bits.length(); i++) { 
			if(bits.charAt(i) == '1') { 
				count += 1; 
			}
		}
		return count;
	}

	public static int calculateU(String bits) { 
		int u = 0;
		for(int i = 0; i < bits.length(); i++) { 
			if(bits.charAt(i) == '1') { 
				for(int j = i + 1; j < bits.length(); j++) { 
					if(bits.charAt(j) == '0') { 
						u++;
					}
				}
			}
		}

		return u;
	}

	public static int calculateUFromT(int m, int n, String bits) {
		int t = 0;
		for(int i = 0; i < bits.length(); i++) { 
			if(bits.charAt(i) == '0') { 
				t += (i+1);
			}
		}

		return (m*n) + (m * (m + 1)) / 2 - t;
	}

	public double getProb(int m, int n, int u) { 
		int[] argArray = new int[3];
		argArray[0] = m; 
		argArray[1] = n; 
		argArray[2] = u;
		return getValue(new Arguments(argArray));
	}

	public double getLowerPValue(int m, int n, int u) { 
		int[] argArray = new int[3];
		argArray[0] = m; argArray[1] = n;
		double sum = 0.0;
		for(int i = 0; i <= u; i++) {
			argArray[2] = i;
			sum += getValue(new Arguments(argArray));
		}
		return sum;
	}

	public double getUpperPValue(int m, int n, int u) { 
		int[] argArray = new int[3];
		argArray[0] = m; argArray[1] = n;
		int max = m * n;
		double sum = 0.0;
		for(int i = u; i <= max; i++) {
			argArray[2] = i;
			sum += getValue(new Arguments(argArray));
		}
		return sum;
	}

	/* (non-Javadoc)
	 * @see edu.mit.csail.psrg.utils.RecurrenceEquation#getValue(edu.mit.csail.psrg.utils.RecurrenceEquation.Arguments)
	 */
	public double getValue(Arguments a) {
		if(cache.containsKey(a)) { return cache.get(a); }

		int n = a.getArgument(0), m = a.getArgument(1), u = a.getArgument(2);

		if(u < 0) { 
			return 0.0;
		}

		if(n == 0 || m == 0) { 
			if(u == 0) { 
				return 1.0;
			} else { 
				return 0.0;
			}
		}

		int[] larray = {n-1, m, u-m};
		int[] rarray = {n, m-1, u};
		Arguments left = new Arguments(larray);
		Arguments right = new Arguments(rarray);
		double lf = (double)n / (double)(n + m);
		double rf = (double)m / (double)(n + m);
		double leftValue = getValue(left);
		double rightValue = getValue(right);

		double val = lf * leftValue + rf * rightValue;
		cache.put(a, val);

		/*
        System.out.print("[" + n + "," + m + "," + u + "]");
        System.out.print(" = ");
        System.out.print("(" + n + "/" + (n+m) + ") * ");
        System.out.print(leftValue);
        System.out.print(" + (" + m + "/" + (n+m) + ") * ");
        System.out.print(rightValue);
        System.out.println(" = " + val);
		 */

		return val;
	}

	public void printCache(PrintStream p) {
		for (Arguments a : cache.keySet()) {
			p.println(a.getArgument(0)+"\t"+a.getArgument(1)+"\t"+a.getArgument(2)+"\t"+
					cache.get(a));
		}
	}

}
