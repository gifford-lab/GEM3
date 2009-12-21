package edu.mit.csail.cgs.clustering.affinitypropagation;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import edu.mit.csail.cgs.clustering.Clusterable;

/**
 * 
 * @author reeder
 *
 */
public class APCluster {

	public static double cluster(Vector<Clusterable> objects, SimilarityMeasure<Clusterable> s, double lam, int convit, int maxit) {

		s.addNoise();
		HashMap<Pair, Double> a = new HashMap<Pair, Double>();
		HashMap<Pair, Double> r = new HashMap<Pair, Double>();
		double max1 = SimilarityMeasure.NEGINF, max2 = SimilarityMeasure.NEGINF;
		int i1=0;
		double tmp = 0.0;
		int it = 0, decit = 0;
		boolean done = false;
		int[][] e = new int[s.size()][convit];
		int[] se = new int[s.size()];
		int decsumc = 0, decsum0 = 0;
		
		//initialize variables
		for (int i=0; i<s.size(); i++) {
			for (int j=0; j<s.size(); j++) {
				Pair ijpair = new Pair(objects.get(i), objects.get(j));
				if (s.exists(ijpair)) {
					a.put(ijpair, 0.0);
					r.put(ijpair, 0.0);
				}
			}
			for (int j=0; j<convit; j++) {
				e[i][j] = 0;
			}
			se[i] = 0;
		}
		
		System.err.println("a size: "+a.size());
		System.err.println("r size: "+r.size());
		
		/*
		System.out.println("a:");
		print(a);
		System.out.println();
		System.out.println("r:");
		print(r);
		System.out.println();
		*/

		while (!done) {
			if (it%10 == 0) {
				System.err.println("Iteration: "+it);
			}
			
			//compute responsibilities
			for (int i=0; i<s.size(); i++) {
				max1 = SimilarityMeasure.NEGINF; max2 = SimilarityMeasure.NEGINF;
				i1=0;
				for (int k=0; k<s.size(); k++) {
					Pair ikpair = new Pair(objects.get(i), objects.get(k));
					if (!s.exists(ikpair)) {
						//System.out.println("!exists");
						continue;
					}
					tmp = a.get(ikpair) + s.evaluate(ikpair);
					//if (it==0) System.out.println(i+" "+k+" "+tmp);
					if (tmp > max1) {
						max2 = max1;
						max1 = tmp;
						i1 = k;
					} else if (tmp > max2) {
						max2 = tmp;
					}
				}
				System.err.println("MAX1: " + max1 + " MAX2: " + max2);
				for (int k=0; k<s.size(); k++) {
					Pair ikpair = new Pair(objects.get(i), objects.get(k));
					if (!s.exists(ikpair)) continue;
					if (k==i1) {
						double message = lam*r.get(ikpair) + (1.0-lam)*(s.evaluate(ikpair) - max2);
						System.err.println("R " + objects.get(i).name() + " to " + objects.get(k).name() + ": " + message);
						r.put(ikpair, message);
					} else {
						double message = lam*r.get(ikpair) + (1.0-lam)*(s.evaluate(ikpair) - max1);
						System.err.println("R " + objects.get(i).name() + " to " + objects.get(k).name() + ": " + message);
						r.put(ikpair, message);
					}
				}
			}
			/*
			if (false) {
				System.out.println("max1: "+max1);
				System.out.println("max2: "+max2);
				System.out.println("i1: "+i1);
				for (int i=0; i<s.size(); i++) {
					for (int j=0; j<s.size(); j++) {
						if (!s.exists(objects.get(i), objects.get(j))) continue;
						System.out.print(r.get(objects.get(i).name()+"SEP"+objects.get(j).name()));
						System.out.print("\t");
					}
					System.out.println();
				}
			}
			
			/*
			if (it==0) {
				System.out.println("max1: "+max1);
				System.out.println("max2: "+max2);
				System.out.println("r("+it+"):");
				print(r);
				System.out.println();
			}
			*/

			//compute availabilities
			for (int k=0; k<s.size(); k++) {
				tmp = 0.0;
				for (int i=0; i<s.size(); i++) {
					Pair ikpair = new Pair(objects.get(i), objects.get(k));
					if (!s.exists(ikpair)) continue;
					if ((i!=k)&&(r.get(ikpair)>0)) {
						tmp += r.get(ikpair);
					}
				}
				for (int i=0; i<s.size(); i++) {
					Pair ikpair = new Pair(objects.get(i), objects.get(k));
					if (!s.exists(ikpair)) continue;
					double tmp2 = tmp;
					double ikr = r.get(ikpair);
					if ((i!=k)&&(ikr>0)) {
						tmp2 -= ikr;
					}
					if (i!=k) {
						tmp2 += r.get(new Pair(objects.get(k),objects.get(k)));
					}
					if (i==k) {
						double message = lam*a.get(ikpair) + (1.0 - lam)*tmp2;
						System.err.println("A " + objects.get(i).name() + " to " + objects.get(k).name() + ": " + message);
						a.put(ikpair, message);
					} else if (tmp2 < 0) {
						double message = lam*a.get(ikpair) + (1.0-lam)*tmp2;
						System.err.println("A " + objects.get(i).name() + " to " + objects.get(k).name() + ": " + message);
						a.put(ikpair, message);
					} else {
						double message = lam*a.get(ikpair);
						System.err.println("A " + objects.get(i).name() + " to " + objects.get(k).name() + ": " + message);
						a.put(ikpair, message);
					}
				}
			}

			//check for convergence
			for (int j = 0; j<s.size(); j++) {
				Pair jjpair = new Pair(objects.get(j), objects.get(j));
				if (a.get(jjpair)+r.get(jjpair) > 0) {
					e[j][decit] = 1;
				} else {
					e[j][decit] = 0;
				}
				se[j] = 0;
				for (int d=0; d<convit; d++) {
					se[j] += e[j][d];
				}
			}
			
			/*
			System.out.println("e("+it+"):");
			print(e);
			*/
			
			decsumc = 0;
			decsum0 = 0;
			for (int j = 0; j<s.size(); j++) {
				if (se[j]==convit) {
					decsumc++;
				} else if (se[j]==0) {
					decsum0++;
				}
			}
			if (((decsumc+decsum0==s.size())&&(decsumc>0))||(it>=maxit)) {
				done = true;
				decit--;
			}
			
			it++; decit++;
			if (decit>=convit) {
				decit = 0;
			}
		}

		System.out.println("iterations: "+it);
		//compute assignments to exemplars
		//e[][decit] represents the exemplars
		//decsumc is the number of exemplars
		Vector<Integer> exidx = new Vector<Integer>();
		int tmpidx = 0;
		int[] assgn = new int[s.size()];
		
		for (int i=0; i<s.size(); i++) {
			if (e[i][decit]==1) {
				exidx.add(i);
			}
		}
		for (int i=0; i<s.size(); i++) {
			max1 = SimilarityMeasure.NEGINF;
			for (int j=0; j<exidx.size(); j++) {
				Pair ijpair = new Pair(objects.get(i), objects.get(exidx.get(j)));
				if (s.evaluate(ijpair) > max1) {
					max1 = s.evaluate(ijpair);
					assgn[i] = j;
				}
			}
		}
		for (int i=0; i<decsumc; i++) {
			assgn[exidx.get(i)] = i;
		}
		
		int[] intarry = new int[1];
		
		s.putAssignments(assgn);
		s.putExemplars(toIntArray(exidx));
		
		double netsim = 0.0;
		for (int i=0; i<s.size(); i++) {
    		netsim += s.evaluate(new Pair(objects.get(i), objects.get(exidx.get(assgn[i]))));
    	}
		
		return netsim;
	}
	
	public static int[] toIntArray(Vector<Integer> vec) {
		int[] toreturn = new int[vec.size()];
		for (int i=0; i<toreturn.length; i++) {
			toreturn[i] = vec.get(i);
		}
		return toreturn;
	}
	
	private static void print(double[][] a) {
		for (int i=0; i<a.length; i++) {
			for (int j=0; j<a[i].length; j++) {
				System.out.print(a[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	private static void print(int[][] a) {
		for (int i=0; i<a.length; i++) {
			for (int j=0; j<a[i].length; j++) {
				System.out.print(a[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	/**
	 * argv[0] = similarities file
	 * argv[1] = preference value
	 * argv[2] = results file
	 * argv[3] = indices file
	 * @param argv
	 * @throws Exception
	 */
	public static void main(String[] argv) throws Exception {
		FileSimilarityMeasure<Clusterable> fsm = new FileSimilarityMeasure<Clusterable>(argv[0],
				"  ", Double.valueOf(argv[1]).doubleValue());
		double netsim = cluster(fsm.objects(), fsm, 0.5, 50, 500);
		PrintStream outstream = new PrintStream(argv[2]);
		fsm.printExemplars(outstream);
		outstream.println();
		fsm.printAssignments(outstream);
		outstream.println();
		outstream.println("Net Similarity: "+netsim);
		fsm.printClusterCenterIndices(new PrintStream(argv[3]));
	}

}