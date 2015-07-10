package edu.mit.csail.cgs.deepseq.discovery.kmer.mtree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;

public class ProbMer {
	
	public double a;
	public double t;
	public double c;
	public double g;
	
	private final ArrayList<Character> monomers = (ArrayList<Character>) Arrays.asList('A', 'T', 'C', 'G');
	
	public ProbMer(double a, double t, double c, double g) {
		this.a = a;
		this.t = t;
		this.c = c;
		this.g = g;
	}
	
	public static double ycDistance(List<ProbMer> kmer1, List<ProbMer> kmer2) {
		if (kmer1.size() == kmer2.size()) {
			double dist = 0;
			for (int i = 0; i < kmer1.size(); i++) {
				dist += (Math.abs(kmer1.get(i).a - kmer2.get(i).a) + Math.abs(kmer1.get(i).t - kmer2.get(i).t)
						+ Math.abs(kmer1.get(i).c - kmer2.get(i).g) + Math.abs(kmer1.get(i).g - kmer2.get(i).g));
			}
			return dist/2;
		}
		else if (kmer1.size() < kmer2.size()) {
			double dist = kmer1.size();
			for (int i = 0; i < kmer2.size() - kmer1.size() + 1; i++) {
				dist = Math.min(dist, ProbMer.ycDistance(kmer1, kmer2.subList(i,  i + kmer1.size())));
			}
			return dist + kmer2.size() - kmer1.size();
		}
		else {
			return ProbMer.ycDistance(kmer2, kmer1);
		}
	}
	
	public static double ycDistance(String kmer1, String kmer2) {
		return CommonUtils.strMinDistance(kmer1, kmer2);
	}
	
	public static ArrayList<ProbMer> rc(ArrayList<ProbMer> forward) {
		ArrayList<ProbMer> rc = new ArrayList<ProbMer>();
		for (int i = 0; i < forward.size(); i++) {
			ProbMer back = forward.get(forward.size() - i - 1);
			rc.add(new ProbMer(1 - back.t, 1 - back.a, 1 - back.g, 1 - back.c));
		}
		return rc;
	}
	
	public static ArrayList<ProbMer> average(ArrayList<ArrayList<ProbMer>> kmerset) {
		ArrayList<ProbMer> average = new ArrayList<ProbMer>();
		int n = kmerset.size();
		for (int i = 0; i < kmerset.get(0).size(); i++) {
			double avga = 0;
			double avgt = 0;
			double avgc = 0;
			double avgg = 0;
			for (List<ProbMer> kmer : kmerset) {
				avga += kmer.get(i).a;
				avgt += kmer.get(i).t;
				avgc += kmer.get(i).c;
				avgg += kmer.get(i).g;
			}
			avga = avga/n;
			avgt = avgt/n;
			avgc = avgc/n;
			avgg = avgg/n;
			average.add(new ProbMer(avga, avgt, avgc, avgg));
		}
		return average;
	}
	
	public ArrayList<ArrayList<ProbMer>> generateKmers(ArrayList<String> kmers, int n) {
		ArrayList<String> kmers0 = new ArrayList<String>(); // no gapmers
		ArrayList<String> kmers1 = new ArrayList<String>(); // 1 gapmer
		ArrayList<String> kmers2 = new ArrayList<String>(); // 2 gapmers
		ArrayList<String> kmers3 = new ArrayList<String>(); // 3 gapmers
		ArrayList<ArrayList<ProbMer>> weighted = new ArrayList<ArrayList<ProbMer>>();
		for (String s: kmers) {
			if (s.length() == n) {
				kmers0.add(s);
			}
			else if (s.length() == n+1) {
				kmers1.add(s);
			}
			else if (s.length() == n+2) {
				kmers2.add(s);
			}
			else {
				kmers3.add(s);
			}
		}
		for (String s: kmers0) {
			ArrayList<ProbMer> weights = new ArrayList<ProbMer>();
			for (int i = 0; i < s.length(); i++) {
				switch (s.charAt(i)) {
					case 'A':
						weights.add(new ProbMer(1, 0, 0, 0));
						break;
					case 'T':
						weights.add(new ProbMer(0, 1, 0, 0));
						break;
					case 'C':
						weights.add(new ProbMer(0, 0, 1, 0));
						break;
					case 'G':
						weights.add(new ProbMer(0, 0, 0, 1));
						break;
					default:
						weights.add(new ProbMer(0.25, 0.25, 0.25, 0.25));
						break;
				}
			}
			weighted.add(weights);
		}
		HashMap<Integer, ProbMer> stored1 = new HashMap<Integer, ProbMer>();
		for (String s: kmers1) {
			ArrayList<ProbMer> weights = new ArrayList<ProbMer>();
			for (int i = 0; i < s.length(); i++) {
				switch (s.charAt(i)) {
					case 'A':
						weights.add(new ProbMer(1, 0, 0, 0));
						break;
					case 'T':
						weights.add(new ProbMer(0, 1, 0, 0));
						break;
					case 'C':
						weights.add(new ProbMer(0, 0, 1, 0));
						break;
					case 'G':
						weights.add(new ProbMer(0, 0, 0, 1));
						break;
					case 'N':
						if (stored1.containsKey(i)) {
							weights.add(stored1.get(i));
						}
						else {
							double acount = 0;
							double tcount = 0;
							double ccount = 0;
							double gcount = 0;
							for (String w: kmers1) {
								switch (w.charAt(i)) {
									case 'A':
										acount++;
										break;
									case 'T':
										tcount++;
										break;
									case 'C':
										ccount++;
										break;
									case 'G':
										gcount++;
										break;
									default:
										break;
								}
							}
							double total = acount + tcount + ccount + gcount;
							if (total != 0) {
								acount /= total;
								tcount /= total;
								ccount /= total;
								gcount /= total;
							}
							else {
								acount = 0.25;
								tcount = 0.25;
								ccount = 0.25;
								gcount = 0.25;
							}
							ProbMer weightedMonomer = new ProbMer(acount, tcount, ccount, gcount);
							weights.add(weightedMonomer);
							stored1.put(i, weightedMonomer);
						}
						break;
					default:
						weights.add(new ProbMer(0.25, 0.25, 0.25, 0.25));
						break;
				}
			}
			weighted.add(weights);
		}
		for (String s: kmers2) {
			ArrayList<ProbMer> weights = new ArrayList<ProbMer>();
			ArrayList<Integer> gapIndices = new ArrayList<Integer>();
			HashMap<ArrayList<Integer>, ArrayList<ProbMer>> stored2 = new HashMap<ArrayList<Integer>, ArrayList<ProbMer>>();
			for (int i = 0; i < s.length(); i++) {
				switch (s.charAt(i)) {
					case 'A':
						weights.add(new ProbMer(1, 0, 0, 0));
						break;
					case 'T':
						weights.add(new ProbMer(0, 1, 0, 0));
						break;
					case 'C':
						weights.add(new ProbMer(0, 0, 1, 0));
						break;
					case 'G':
						weights.add(new ProbMer(0, 0, 0, 1));
						break;
					case 'N':
						weights.add(new ProbMer(0.25, 0.25, 0.25, 0.25));
						gapIndices.add(i);
						break;
					default:
						weights.add(new ProbMer(0.25, 0.25, 0.25, 0.25));
						break;
				}
			}
			if (stored2.containsKey(gapIndices)) {
				ArrayList<ProbMer> indexedProbMers = stored2.get(gapIndices);
				for (int j = 0; j < gapIndices.size(); j++) {
					weights.set(gapIndices.get(j), indexedProbMers.get(j));
				}
			}
			else {
				ArrayList<ProbMer> indexedProbMers = new ArrayList<ProbMer>();
				double[][] counts = new double[4][4];
				for (String w: kmers2) {
					if (w.charAt(gapIndices.get(0)) != 'N' && w.charAt(gapIndices.get(1)) != 'N') {
						counts[monomers.indexOf(w.charAt(gapIndices.get(0)))][monomers.indexOf(w.charAt(gapIndices.get(1)))]++;
					}
				}
				
				stored2.put(gapIndices, indexedProbMers);
			}
		}
		return weighted;
	}
	
	
}
