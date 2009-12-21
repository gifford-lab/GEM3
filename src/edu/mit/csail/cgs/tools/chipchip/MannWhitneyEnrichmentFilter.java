package edu.mit.csail.cgs.tools.chipchip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

import com.hp.hpl.jena.util.iterator.ConcatenatedIterator;
import com.hp.hpl.jena.util.iterator.NullIterator;

import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.Filter;
import edu.mit.csail.cgs.ewok.verbs.FilterIterator;
import edu.mit.csail.cgs.ewok.verbs.GeneToPromoter;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.SetTools;
import edu.mit.csail.cgs.utils.probability.MannWhitneyEquation;
import edu.mit.csail.cgs.utils.probability.Sample;
import edu.mit.csail.cgs.datasets.species.Gene;

public class MannWhitneyEnrichmentFilter<X extends Region> implements Filter<X,X> {

	Genome genome;
	public String mark, stage, version;
	private ChipChipLocator loc;
	private ChipChipData data;
	private MannWhitneyEquation eq;
	private double pValueCutoff;
	private static final String[] geneTables = {"ensGene", "mgcGenes", "refGene"};
	private static final boolean PCORRECTION = true;
	private RandomEngine engine;
	private Normal norm;

	public MannWhitneyEnrichmentFilter(Genome g, String m, String s, String v,
			double pValueCutoff, MannWhitneyEquation eq) {
		genome = g;
		mark = m; stage = s;
		String name;
		if (m.equals("H3K79me2") && (s.equals("ES+2d Stage") || s.equals("ES+2d Stage, 8 hours post RA") || s.equals("2+1 day"))) {
			name = "Mm " + mark + ":HBG3:" + stage + 
			" vs WCE:HBG3:" + stage;
		} else {
			name = "Mm " + mark + ":HBG3:" + stage + 
			" vs H3:HBG3:" + stage;
		}
		version = v;

		loc = new ChipChipLocator(genome, name, version);
		data = loc.createObject();
		this.eq = eq;
		this.pValueCutoff = pValueCutoff;
		engine = new cern.jet.random.engine.DRand(); 
		norm = new Normal(0.0, 1.0, engine);
	}
	
	public boolean equals(Object o) {
		if (o instanceof MannWhitneyEnrichmentFilter) {
			MannWhitneyEnrichmentFilter other = (MannWhitneyEnrichmentFilter)o;
			return this.mark.equals(other.mark) && this.stage.equals(other.stage) && this.version.equals(other.version);
		} else {
			return false;
		}
	}

	/*
	public MannWhitneyEnrichmentFilter(Genome g, String m, String s, String v,
			double pValueCutoff, File cacheFile) {
		genome = g;
		mark = m; stage = s;
		String name = "Mm " + mark + ":HBG3:" + stage + 
		" vs H3:HBG3:" + stage;
		String version = v;

		loc = new ChipChipLocator(genome, name, version);
		data = loc.createObject();
		eq = new MannWhitneyEquation(cacheFile);
		this.pValueCutoff = pValueCutoff;
	}
	 */

	/**
	 * args[0] = genome
	 * args[1] = mark
	 * args[2] = stage
	 * args[3] = version
	 * args[4] = p-value cutoff
	 * args[5] = region file
	 * args[6] = output file
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		/*
		Genome mm8 = Organism.findGenome("mm8");
		ComparableGene g1 = new ComparableGene(mm8, "1", 16086510, 16089593, "BC086786", "1", "ref");
		ComparableGene g2 = new ComparableGene(mm8, "1", 16086510, 16089597, "BC025909", "1", "ref");
		Set<ComparableGene> testset = new HashSet<ComparableGene>();
		testset.add(g1);
		System.err.println(g1.equals(g2));
		 */
		/*
		Genome mm8 = Organism.findGenome("mm8");
		String mark = "H3K27me3";
		String stage = "ES Stage";
		String version = "median linefit";
		double pValueCutoff = Double.parseDouble("0.05");

		List<NamedStrandedRegion> promList = promListFromFile(new File(args[5]), mm8);

		if (PCORRECTION) {
			pValueCutoff /= ((double)promList.size());
		}
		MannWhitneyEquation eq = new MannWhitneyEquation();
		Iterator<NamedStrandedRegion> promIter = new FilterIterator<NamedStrandedRegion,NamedStrandedRegion>(new MannWhitneyEnrichmentFilter<NamedStrandedRegion>(mm8,mark,stage,version,pValueCutoff, eq),
				promList.iterator());


		PrintStream out = new PrintStream(args[6]);
		System.err.println("writing enriched gene file");
		while (promIter.hasNext()) {
			NamedStrandedRegion tmpProm = promIter.next();
			if (tmpProm!=null) {
				out.println(tmpProm.getName()+"\t"+tmpProm.getChrom()+"\t"+tmpProm.getStart()+
						"\t"+tmpProm.getEnd()+"\t"+tmpProm.getStrand());
			}
		}
		out.flush();
		out.close();
		 */

		Genome mm8 = Organism.findGenome("mm8");
		String mark = "H3K27me3";
		String stage = "ES+2d Stage, before RA";
		String version = "median linefit";
		MannWhitneyEquation eq = new MannWhitneyEquation();
		double pValueCutoff = Double.parseDouble("0.05");
		PrintStream out = new PrintStream(args[0]);
		MannWhitneyEnrichmentFilter<NamedStrandedRegion> mwef = new MannWhitneyEnrichmentFilter<NamedStrandedRegion>(mm8, mark, stage, version,
				pValueCutoff, eq);
		mwef.printMethodComparison(50, 50, 1000, out);
		out.flush();
		out.close();
	}

	public static Iterator<NamedStrandedRegion> promIterFromFile(File f, Genome g) throws IOException {
		LinkedList<NamedStrandedRegion> promList = new LinkedList<NamedStrandedRegion>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		String[] split;
		while((line = br.readLine()) != null) { 
			split = line.split("\t");
			promList.addLast(new NamedStrandedRegion(g, split[1], Integer.parseInt(split[2]), Integer.parseInt(split[3]), split[0], split[4].charAt(0)));
		}
		br.close();
		return promList.iterator();
	}

	public static List<NamedStrandedRegion> promListFromFile(File f, Genome g) throws IOException {
		LinkedList<NamedStrandedRegion> promList = new LinkedList<NamedStrandedRegion>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		String[] split;
		while((line = br.readLine()) != null) { 
			split = line.split("\t");
			promList.addLast(new NamedStrandedRegion(g, split[1], Integer.parseInt(split[2]), Integer.parseInt(split[3]), split[0], split[4].charAt(0)));
		}
		br.close();
		return promList;
	}
	
	public static List<NamedStrandedRegion> promListFromTSSFile(File f, Genome g, int upstream, int downstream) throws IOException {
		LinkedList<NamedStrandedRegion> promList = new LinkedList<NamedStrandedRegion>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		String[] split;
		while((line = br.readLine()) != null) { 
			split = line.split("\t");
			int start, stop;
			switch(split[3].charAt(0)) { 
	        default:
	        case '+':
	            start = Integer.parseInt(split[2]) - upstream;
	            stop = Integer.parseInt(split[2]) + downstream;
	        case '-':
	            start = Integer.parseInt(split[2]) - downstream;
	            stop = Integer.parseInt(split[2]) + upstream;
	        }
			promList.addLast(new NamedStrandedRegion(g, split[1], start, stop, split[0], split[3].charAt(0)));
		}
		br.close();
		return promList;
		
        
	}

	public X execute(X a) {
		double pvalue = computePValue(a);
		if (pvalue <= pValueCutoff) {
			return a;
		} else {
			return null;
		}
	}

	public double computePValue(X a) {
		try {
			data.window(a.getChrom(), a.getStart(), a.getEnd());
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<Sample> wceSamples = new ArrayList<Sample>();
		ArrayList<Sample> ipSamples = new ArrayList<Sample>();
		for (int i=0; i<data.getCount(); i++) {
			for (int j=0; j<data.getReplicates(i); j++) {
				wceSamples.add(new Sample(data.getWCE(i, j),0));
				ipSamples.add(new Sample(data.getIP(i, j),1));
			}
		}
		if (true) {
			int u = computeU(wceSamples, ipSamples);
			System.err.println(ipSamples.size()+"\t"+wceSamples.size()+"\t"+u);
			return eq.getLowerPValue(ipSamples.size(), wceSamples.size(), u);
		} else {
			return computeNormalApprox(wceSamples, ipSamples);
		}
	}

	private void printMethodComparison(int mMax, int nMax, int reps, PrintStream out) {
		Random r = new Random();
		for (int m=1; m<=mMax; m++) {
			for (int n=1; n<=nMax; n++) {
				double[] diffs = new double[reps];
				for (int rep=0; rep<reps; rep++) {
					diffs[rep] = Math.abs(compareMethods(m, n, r));
				}
				out.println(m+"\t"+n+"\t"+average(diffs));
			}
		}
	}

	public static double average(double[] arr) {
		double sum = 0.0d;
		for (int i=0; i<arr.length; i++) {
			sum += arr[i];
		}
		return sum / ((double)arr.length);
	}

	private double compareMethods(int m, int n, Random r) {
		List<Sample> controlSamples = new ArrayList<Sample>();
		List<Sample> testSamples = new ArrayList<Sample>();
		double value = 0.0d;
		double tmp = r.nextDouble();
		int mIndex = m;
		int nIndex = n;
		while (mIndex>0 || nIndex>0) {
			if (r.nextDouble()>tmp) {
				if (mIndex>0) {
					testSamples.add(new Sample(value, 1));
					mIndex--;
				} else {
					controlSamples.add(new Sample(value, 0));
					nIndex--;
				}
			} else {
				if (nIndex>0) {
					controlSamples.add(new Sample(value, 0));
					nIndex--;
				} else {
					testSamples.add(new Sample(value, 1));
					mIndex--;
				}
			}
			value += 1.0d;
		}
		int u = computeU(controlSamples, testSamples);
		double normApprox = computeNormalApprox(controlSamples, testSamples);
		double exactVal = eq.getLowerPValue(m, n, u);
		System.err.println(m+"\t"+n+"\t"+normApprox+"\t"+exactVal);
		return normApprox - exactVal;
	}

	private double computeNormalApprox(List<Sample> controlSamples, List<Sample> testSamples) {
		double m = ((double)testSamples.size());
		double n = ((double)controlSamples.size());
		double E = m*(m+n+1.0d)/2.0d;
		double sigma = Math.sqrt(m*n*(m+n+1.0d)/12.0d);
		double t = (double)computeT(controlSamples, testSamples);
		return 1.0d - norm.cdf((t-E)/sigma);
	}

	private static int computeT(List<Sample> controlSamples, List<Sample> testSamples) {
		List<Sample> allSamples = new ArrayList<Sample>(controlSamples);
		allSamples.addAll(testSamples);
		Collections.sort(allSamples);
		int rank = 1;
		int t = 0;
		for (Sample s : allSamples) {
			if (s.getGroup() == 1) {
				t += rank;
			}
			rank++;
		}
		return t;
	}

	private static int computeU(List<Sample> controlSamples, List<Sample> testSamples) {
		int m = testSamples.size();
		int n = controlSamples.size();
		int u = m*n + m*(m+1)/2;
		int t = computeT(controlSamples, testSamples);
		u -= t;
		//System.err.println(m+"\t"+n+"\t"+u);
		return u;
	}

	public void printCache(PrintStream p) {
		eq.printCache(p);
	}



}

class WellTiledGeneGenerator implements Iterator<Gene> { 

	private final String[] geneTables = {"ensGene", "mgcGenes", "refGene"};
	private Genome genome;
	private File file;

	private LinkedList<Gene> tiledGenes;
	private Mapper<Gene,NamedStrandedRegion> gtp;

	public WellTiledGeneGenerator(File file, int tssUp, int tssDown) { 
		try { 
			genome = Organism.findGenome("mm8");
		} catch(NotFoundException nfe) { 
			throw new IllegalStateException();
		}
		this.file = file;
		this.tiledGenes = new LinkedList<Gene>();
		gtp = new GeneToPromoter(tssUp,tssDown);
		parse();
	}

	public void parse() {
		SortedSet<ComparableGene> geneSet = new TreeSet<ComparableGene>();
		for (String table : geneTables) {
			Expander<NamedRegion,Gene> gen = new RefGeneGenerator<NamedRegion>(genome, table);
			Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
			Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
			while (genes.hasNext()) {
				geneSet.add(new ComparableGene(genes.next()));
			}
		}
		System.err.println("number of unique genes: "+geneSet.size());
		SortedSet<Region> regionSet = new TreeSet<Region>();
		try { 
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = null;
			String[] split;
			while((line = br.readLine()) != null) { 
				split = line.split("\t");
				String chrom = split[0];
				int start = Integer.parseInt(split[1]);
				int end = Integer.parseInt(split[2]);
				Region r = new Region(genome, chrom, start, end);
				regionSet.add(r);
			}
			br.close();
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
		Iterator<ComparableGene> geneIter = geneSet.iterator();
		Iterator<Region> regionIter = regionSet.iterator();
		Region tmpRegion = new Region(genome,"",-2,-1);
		while (geneIter.hasNext()) {
			Gene tmpGene = geneIter.next();
			NamedStrandedRegion prom = gtp.execute(tmpGene);
			while (regionIter.hasNext() && (tmpRegion.getEnd()<prom.getEnd() || 
					!tmpRegion.getChrom().equals(prom.getChrom())))
				tmpRegion = regionIter.next();
			if (tmpRegion.contains(prom)) {
				tiledGenes.addLast(tmpGene);
				System.err.println(tmpGene.getName()+"\t"+tmpGene.getChrom()+"\t"+tmpGene.getStart()+
						"\t"+tmpGene.getEnd()+"\t"+tmpGene.getStrand());
			}
		}
	}

	public int getCurrentSize() { return tiledGenes.size(); }
	public boolean hasNext() { return !tiledGenes.isEmpty(); }
	public Gene next() { return tiledGenes.removeFirst(); }
	public void remove() { throw new UnsupportedOperationException(); }
}

class ComparableGene extends Gene {
	public ComparableGene(Genome g, String c, int start, int end, String name, String id, String src) {
		super(g,c,start,end,name,id,src);
	}

	public ComparableGene(Genome g, String c, int start, int end, String name, String id, char str, String src) {
		super(g,c,start,end,name,id,str,src);
	}

	public ComparableGene(Genome g, String c, int start, int end, String name, String id, char str, String src, int dbid) {
		super(g,c,start,end,name,id,str,src,dbid);
	}

	public ComparableGene(Gene g) {
		super(g.getGenome(),g.getChrom(),g.getStart(),g.getEnd(),g.getName(),g.getID(),g.getStrand(),g.getSource(),g.getDBID());
	}

	public boolean equals(Object o) {
		if(!(o instanceof ComparableGene)) { return false; }
		ComparableGene g = (ComparableGene)o;
		return this.getChrom().equals(g.getChrom()) &&
		((this.getStart() == g.getStart()) ||
				(this.getEnd() == g.getEnd()));
	}
}

class CoveredFilter<X extends Region> implements Filter<X,X> {
	private LinkedList<Region> regions;
	private File file;
	private Genome genome;

	public CoveredFilter(File file) {
		try { 
			genome = Organism.findGenome("mm8");
		} catch(NotFoundException nfe) { 
			throw new IllegalStateException();
		}
		this.file = file;
		this.regions = new LinkedList<Region>();
		parse();
	}

	private void parse() {
		try { 
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = null;
			String[] split;
			while((line = br.readLine()) != null) { 
				split = line.split("\t");
				String chrom = split[0];
				int start = Integer.parseInt(split[1]);
				int end = Integer.parseInt(split[2]);
				Region r = new Region(genome, chrom, start, end);
				regions.add(r);
			}
			br.close();
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}

	public X execute(X a) {
		for (Region r : regions) {
			if (r.contains(a)) {
				return a;
			}
		}
		return null;
	}
}

class UniqueishRegionFilter<X extends Region> implements Filter<X,X> {

	private int buffer;
	private Map<List<String>,Set<X>> leftSeen;
	private Map<List<String>,Set<X>> rightSeen;
	private SetTools<X> st;

	public UniqueishRegionFilter(int buffer) {
		this.buffer = buffer;
		this.leftSeen = new HashMap<List<String>,Set<X>>();
		this.rightSeen = new HashMap<List<String>,Set<X>>();
		this.st = new SetTools<X>();
	}

	public X execute(X a) {
		List<String> tmpList = new ArrayList<String>();
		tmpList.add(a.getChrom());
		Set<X> leftSet = new HashSet<X>();
		for (int l=a.getStart()-buffer; l<=a.getStart()+buffer; l++) {
			tmpList.add((new Integer(l)).toString());
			if (leftSeen.containsKey(tmpList)) {
				leftSet.addAll(leftSeen.get(tmpList));
			}
			tmpList.remove(1);
		}
		Set<X> rightSet = new HashSet<X>();
		for (int r=a.getEnd()-buffer; r<=a.getEnd()+buffer; r++) {
			tmpList.add((new Integer(r)).toString());
			if (rightSeen.containsKey(tmpList)) {
				rightSet.addAll(rightSeen.get(tmpList));
			}
			tmpList.remove(1);
		}
		if (st.intersects(leftSet,rightSet)) {
			return null;
		}
		tmpList.add((new Integer(a.getStart())).toString());
		if (!leftSeen.containsKey(tmpList)) {
			leftSeen.put(tmpList, new HashSet<X>());
		}
		leftSeen.get(tmpList).add(a);
		tmpList.remove(1);
		tmpList.add((new Integer(a.getEnd())).toString());
		if (!rightSeen.containsKey(tmpList)) {
			rightSeen.put(tmpList, new HashSet<X>());
		}
		rightSeen.get(tmpList).add(a);
		return a;
	}

}
