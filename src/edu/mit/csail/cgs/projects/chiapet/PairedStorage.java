package edu.mit.csail.cgs.projects.chiapet;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.projects.readdb.Client;
import edu.mit.csail.cgs.projects.readdb.ClientException;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class PairedStorage {

	private Genome g;
	private Map<Integer,String> revChromMap;
	private Map<PairedHit,Integer> index;
	private Map<PairedHit,Double> fullMap;
	private Map<PairedHit,Double> restrictedMap;
	private SortedSet<PairedHit> leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

		public int compare(PairedHit arg0, PairedHit arg1) {
			int tor;
			if (arg0.leftChrom==arg1.leftChrom) {
				tor = arg0.leftPos-arg1.leftPos;
			} else {
				tor = arg0.leftChrom - arg1.leftChrom;
			}
			if (tor==0) {
				if (arg0.rightChrom == arg1.rightChrom) {
					return arg0.rightPos - arg1.rightPos;
				} else {
					return arg0.rightChrom - arg1.rightChrom;
				}
			} else {
				return tor;
			}
		}

	});
	private SortedSet<PairedHit> rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

		public int compare(PairedHit arg0, PairedHit arg1) {
			int tor;
			if (arg0.rightChrom == arg1.rightChrom) {
				tor = arg0.rightPos - arg1.rightPos;
			} else {
				tor = arg0.rightChrom - arg1.rightChrom;
			}
			if (tor==0) {
				if (arg0.leftChrom==arg1.leftChrom) {
					return arg0.leftPos-arg1.leftPos;
				} else {
					return arg0.leftChrom - arg1.leftChrom;
				}
			} else {
				return tor;
			}
		}

	});
	private List<Double> kernel, otherkernel, lik;
	private int halfkern;
	private int minsize, maxsize;
	private int posTotal, negTotal;
	private int chimericReads;
	private Poisson poisson = new Poisson(1, new DRand());
	private int binsize = 100;
	private double tmpnorm = 0.0d;
	private boolean removeArtifact;
	private String dir;
	private String currentChrom = "";

	public static void main(String[] args) throws NotFoundException, IOException, ClientException {
		Genome g = Args.parseGenome(args).cdr();
		String outfile = Args.parseString(args, "outfile", "");
		Region r1 = Region.fromString(g, Args.parseString(args, "region1", ""));
		Region r2 = Region.fromString(g, Args.parseString(args, "region2", ""));
		List<Double> kernel = readDoubleList(Args.parseString(args, "kernel", ""));
		List<Double> otherkernel = readDoubleList(Args.parseString(args, "otherkernel", ""));
		PairedStorage storage = new PairedStorage(g, 0, 0, null, kernel, otherkernel, 0, 0);
		storage.lik = readDoubleList(Args.parseString(args, "lik", ""));
		String hitfile = Args.parseString(args, "hitfile", "");
		String zfile = Args.parseString(args, "zfile", "");
		storage.initializeFromFile(hitfile);
		double[][] z = new double[3][storage.index.size()];
		BufferedReader r = new BufferedReader(new FileReader(zfile));
		String s;
		String[] split;
		int i = 0;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			for (int j=0; j<3; j++) {
				z[j][i] = Double.parseDouble(split[j]);
			}
			i++;
		}
		for (int c=0; c<3; c++) {
			double[][] tor;
			if (Args.parseFlags(args).contains("balanced")) {
				tor = storage.balancedHeatmap(r1, r2, z[c]);
			} else {
				tor = storage.heatmap(r1, r2, z[c]);
			}
			PrintStream out = new PrintStream(outfile+c+".txt");
			for (int c1=0; c1<tor.length; c1++) {
				for (int c2=0; c2<tor[c1].length; c2++) {
					out.print(tor[c1][c2]+"\t");
				}
				out.println();
			}
			out.flush();
			out.close();
		}

	}

	public static List<Double> readDoubleList(String file) throws IOException {
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		List<Double> tor = new ArrayList<Double>();
		while ((s = r.readLine()) != null) {
			tor.add(Double.parseDouble(s));
		}
		return tor;
	}

	public PairedStorage(Genome g, int totalReads, int chimericReads, Set<String> aligns, List<Double> kernel, List<Double> otherkernel, int minsize, int maxsize) {
		this.g = g;
		this.chimericReads = chimericReads;
		this.kernel = kernel;
		this.otherkernel = otherkernel;
		this.halfkern = kernel.size()/2;
		this.minsize = minsize;
		this.maxsize = maxsize;
		/*
		Client client = new Client();
		ChromosomeGenerator<Genome> cg = new ChromosomeGenerator<Genome>();
		Iterator<Region> chromiter = cg.execute(g);
		while (chromiter.hasNext()) {
			Region chrom = chromiter.next();
			try {
				for (String a : aligns) {
					List<PairedHit> hits = client.getPairedHits(a, g.getChromID(chrom.getChrom()), true, chrom.getStart(), chrom.getEnd(), null, null);
					leftset.addAll(hits);
					rightset.addAll(hits);
				}
			} catch (Exception e) {
				System.err.println(e);
			}
		}
		 */
	}
	
	public void removeArtifact() {
		removeArtifact = true;
	}

	public int getIndex(PairedHit hit) {
		return index.get(hit);
	}

	public Set<PairedHit> getHitSet() {
		return index.keySet();
	}

	public Set<PairedHit> getHitSet(Region r) {
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), r.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		Set<PairedHit> subset = new HashSet<PairedHit>();
		subset.addAll(leftset.subSet(from, to));
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), r.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0f);
		subset.addAll(rightset.subSet(from, to));
		return subset;
	}

	public void initializeFromSet(Set<PairedHit> set) {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		for (PairedHit hit : set) {
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
				hit.flipSides();
			}
			if (hit.leftStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			if (hit.rightStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			leftset.add(hit);
			rightset.add(hit);
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}
	
	public boolean isArtifact(PairedHit hit) {
		return hit.lesserStrand() && !hit.greaterStrand() && Math.abs(hit.leftPos-hit.rightPos)<=70;
	}

	public void initializeFromFile(String file) throws NumberFormatException, IOException {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			PairedHit hit = new PairedHit(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Boolean.parseBoolean(split[2]), 
					Short.parseShort(split[3]), Integer.parseInt(split[4]), Integer.parseInt(split[5]), Boolean.parseBoolean(split[6]), 
					Short.parseShort(split[7]), Float.parseFloat(split[8]));
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
				hit.flipSides();
			}
			if (hit.leftStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			if (hit.rightStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			leftset.add(hit);
			rightset.add(hit);
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}
	
	public void initializeFromAnnotatedDirectory(String dir) {
		this.dir = dir;
	}
	
	public String getCurrentChrom() {
		return currentChrom;
	}
	
	public void loadChrom(String chrom) throws NumberFormatException, IOException {
		if (!chrom.equals(currentChrom)) {
			currentChrom = chrom;
			System.err.println("loading "+chrom);
			initializeFromAnnotatedFile(dir+chrom+".txt");
		}
	}
	
	public void initializeFromAnnotatedFile(String file) throws NumberFormatException, IOException {
		fullMap = new HashMap<PairedHit,Double>();
		restrictedMap = new HashMap<PairedHit,Double>();
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			PairedHit hit = new PairedHit(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Boolean.parseBoolean(split[2]), 
					Short.parseShort(split[3]), Integer.parseInt(split[4]), Integer.parseInt(split[5]), Boolean.parseBoolean(split[6]), 
					Short.parseShort(split[7]), Float.parseFloat(split[8]));
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
				hit.flipSides();
			}
			if (hit.leftStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			if (hit.rightStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			leftset.add(hit);
			rightset.add(hit);
			fullMap.put(hit, Double.parseDouble(split[9]));
			restrictedMap.put(hit, Double.parseDouble(split[10]));
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}
	
	public void initializeFromPointFile(String file) throws NumberFormatException, IOException {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			Point left = Point.fromString(g, split[0]);
			Point right = Point.fromString(g, split[1]);
			PairedHit hit = new PairedHit(g.getChromID(left.getChrom()), left.getLocation(), true, 
					Short.valueOf("20"), g.getChromID(right.getChrom()), right.getLocation(), true, 
					Short.valueOf("20"), 1f);
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
				hit.flipSides();
			}
			if (hit.leftStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			if (hit.rightStrand) {
				posTotal++;
			} else {
				negTotal++;
			}
			leftset.add(hit);
			rightset.add(hit);
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}

	public void initializeFromFile(String file, int subclass) throws NumberFormatException, IOException {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			if (Integer.parseInt(split[9])==subclass) {
				PairedHit hit = new PairedHit(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Boolean.parseBoolean(split[2]), 
						Short.parseShort(split[3]), Integer.parseInt(split[4]), Integer.parseInt(split[5]), Boolean.parseBoolean(split[6]), 
						Short.parseShort(split[7]), Float.parseFloat(split[8]));
				if (removeArtifact && isArtifact(hit)) {
					continue;
				}
				if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
					hit.flipSides();
				}
				if (hit.leftStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				if (hit.rightStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				leftset.add(hit);
				rightset.add(hit);
			}
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}

	public void initializeFromFileAnd(String file, Region r1, Region r2) throws NumberFormatException, IOException {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			PairedHit hit = new PairedHit(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Boolean.parseBoolean(split[2]), 
					Short.parseShort(split[3]), Integer.parseInt(split[4]), Integer.parseInt(split[5]), Boolean.parseBoolean(split[6]), 
					Short.parseShort(split[7]), Float.parseFloat(split[8]));
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if ((containsLeft(r1,hit) && containsRight(r2,hit)) || (containsLeft(r2,hit) && containsRight(r1,hit))) {
				if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
					hit.flipSides();
				}
				if (hit.leftStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				if (hit.rightStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				leftset.add(hit);
				rightset.add(hit);
			}
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}

	public void initializeFromFileOr(String file, Region r1, Region r2) throws NumberFormatException, IOException {
		leftset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.leftChrom==arg1.leftChrom) {
					tor = arg0.leftPos-arg1.leftPos;
				} else {
					tor = arg0.leftChrom - arg1.leftChrom;
				}
				if (tor==0) {
					if (arg0.rightChrom == arg1.rightChrom) {
						return arg0.rightPos - arg1.rightPos;
					} else {
						return arg0.rightChrom - arg1.rightChrom;
					}
				} else {
					return tor;
				}
			}

		});
		rightset = new TreeSet<PairedHit>(new Comparator<PairedHit>() {

			public int compare(PairedHit arg0, PairedHit arg1) {
				int tor;
				if (arg0.rightChrom == arg1.rightChrom) {
					tor = arg0.rightPos - arg1.rightPos;
				} else {
					tor = arg0.rightChrom - arg1.rightChrom;
				}
				if (tor==0) {
					if (arg0.leftChrom==arg1.leftChrom) {
						return arg0.leftPos-arg1.leftPos;
					} else {
						return arg0.leftChrom - arg1.leftChrom;
					}
				} else {
					return tor;
				}
			}

		});
		index = new HashMap<PairedHit,Integer>();
		int i = 0;
		BufferedReader r = new BufferedReader(new FileReader(file));
		String s;
		String[] split;
		while ((s = r.readLine()) != null) {
			split = s.split("\t");
			PairedHit hit = new PairedHit(Integer.parseInt(split[0]), Integer.parseInt(split[1]), Boolean.parseBoolean(split[2]), 
					Short.parseShort(split[3]), Integer.parseInt(split[4]), Integer.parseInt(split[5]), Boolean.parseBoolean(split[6]), 
					Short.parseShort(split[7]), Float.parseFloat(split[8]));
			if (removeArtifact && isArtifact(hit)) {
				continue;
			}
			if (containsLeft(r1,hit) || containsRight(r2,hit) || containsLeft(r2,hit) || containsRight(r1,hit)) {
				if (hit.leftPointRegion(g).compareTo(hit.rightPointRegion(g))>0) {
					hit.flipSides();
				}
				if (hit.leftStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				if (hit.rightStrand) {
					posTotal++;
				} else {
					negTotal++;
				}
				leftset.add(hit);
				rightset.add(hit);
			}
		}
		for (PairedHit hit : leftset) {
			index.put(hit, i++);
		}
	}

	public static boolean containsLeft(Region r, PairedHit hit) {
		Genome g = r.getGenome();
		return g.getChromID(r.getChrom())==hit.leftChrom && r.getStart()<=hit.leftPos && r.getEnd()>=hit.leftPos;
	}

	public static boolean containsRight(Region r, PairedHit hit) {
		Genome g = r.getGenome();
		return g.getChromID(r.getChrom())==hit.rightChrom && r.getStart()<=hit.rightPos && r.getEnd()>=hit.rightPos;
	}

	public void scatter(Region r1, Region r2, String file) throws FileNotFoundException {
		System.err.println(leftset.size());
		System.err.println(g.getChromID(r2.getChrom()));
		PrintStream out = new PrintStream(file);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), r2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), r2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		Set<PairedHit> goodset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (r1.contains(p.rightPointRegion(g))) {
				int relright = p.rightPos-r1.getStart();
				int relleft = p.leftPos-r2.getStart();
				out.println(relright+"\t"+relleft);
				goodset.add(p);
			}

		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), r2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), r2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (r1.contains(p.leftPointRegion(g)) && !goodset.contains(p)) {
				int relright = p.rightPos-r2.getStart();
				int relleft = p.leftPos-r1.getStart();
				out.println(relleft+"\t"+relright);
			}
		}
		out.flush();
		out.close();
	}

	public double[][] balancedHeatmap(Region r1, Region r2, double[] z) {
		double[][] tor = new double[r1.getWidth()/binsize+1][r2.getWidth()/binsize+1];
		Region e1 = r1.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		Region e2 = r2.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		Set<PairedHit> goodset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (e1.contains(p.rightPointRegion(g))) {
				int pindex = index.get(p);
				int relright = p.rightPos-r1.getStart();
				int relleft = p.leftPos-r2.getStart();

				for (int i=(Math.max(r1.getStart(), p.rightPos-(binsize*lik.size()/2))-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.rightPos+lik.size()*binsize/2)-r1.getStart())/binsize-1; i++) {
					for (int j=(Math.max(r2.getStart(), p.leftPos-(binsize*lik.size()/2))-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.leftPos+lik.size()*binsize/2)-r2.getStart())/binsize-1; j++) {
						//System.err.println(z[pindex]+"\t"+lik.get(i-(relright/binsize))+"\t"+lik.get(j-(relleft/binsize)));
						tor[i][j] += z[pindex]*lik.get(i-(relright/binsize)+lik.size()/2+1)*lik.get(j-(relleft/binsize)+lik.size()/2+1);
					}
				}

				goodset.add(p);
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			if (e1.contains(p.leftPointRegion(g)) && !goodset.contains(p)) {
				int pindex = index.get(p);
				int relright = p.rightPos-r2.getStart();
				int relleft = p.leftPos-r1.getStart();

				for (int i=(Math.max(r1.getStart(), p.leftPos-(binsize*lik.size()/2))-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.leftPos+lik.size()*binsize/2)-r1.getStart())/binsize-1; i++) {
					for (int j=(Math.max(r2.getStart(), p.rightPos-(binsize*lik.size()/2))-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.rightPos+lik.size()*binsize/2)-r2.getStart())/binsize-1; j++) {
						tor[i][j] += z[pindex]*lik.get(i-(relleft/binsize)+lik.size()/2+1)*lik.get(j-(relright/binsize)+lik.size()/2+1);
					}
				}

			}
		}
		return tor;
	}

	public double[][] heatmap(Region r1, Region r2, double[] z) {
		double[][] tor = new double[r1.getWidth()/binsize+1][r2.getWidth()/binsize+1];
		Region e1 = r1.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		Region e2 = r2.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		Set<PairedHit> goodset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (e1.contains(p.rightPointRegion(g))) {
				int pindex = index.get(p);
				int relright = p.rightPos-r1.getStart();
				int relleft = p.leftPos-r2.getStart();
				if (!p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.rightPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.leftPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							//System.err.println(z[pindex]+"\t"+lik.get(i-(relright/binsize))+"\t"+lik.get(j-(relleft/binsize)));
							tor[i][j] += z[pindex]*lik.get(i-(relright/binsize))*lik.get(j-(relleft/binsize));
						}
					}
				}
				if (!p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.rightPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.leftPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relright/binsize)-i)*lik.get(j-(relleft/binsize));
						}
					}
				}
				if (p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.rightPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.leftPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relright/binsize))*lik.get((relleft/binsize)-j);
						}
					}
				}
				if (p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.rightPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.leftPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relright/binsize)-i)*lik.get((relleft/binsize)-j);
						}
					}
				}
				goodset.add(p);
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			if (e1.contains(p.leftPointRegion(g)) && !goodset.contains(p)) {
				int pindex = index.get(p);
				int relright = p.rightPos-r2.getStart();
				int relleft = p.leftPos-r1.getStart();
				if (!p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.leftPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.rightPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relleft/binsize))*lik.get(j-(relright/binsize));
						}
					}
				}
				if (!p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.leftPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.rightPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relleft/binsize)-i)*lik.get(j-(relright/binsize));
						}
					}
				}
				if (p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.leftPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.rightPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relleft/binsize))*lik.get((relright/binsize)-j);
						}
					}
				}
				if (p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.leftPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.rightPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relleft/binsize)-i)*lik.get((relright/binsize)-j);
						}
					}
				}
			}
		}
		return tor;
	}

	public double[][] probabilityMap(Region r1, Region r2, double[] z) {
		double[][] tor = new double[r1.getWidth()/binsize+1][r2.getWidth()/binsize+1];
		Region e1 = r1.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		Region e2 = r2.expand((lik.size()-1)*binsize, (lik.size()-1)*binsize);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		Set<PairedHit> goodset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (e1.contains(p.rightPointRegion(g))) {
				int pindex = index.get(p);
				int relright = p.rightPos-r1.getStart();
				int relleft = p.leftPos-r2.getStart();
				if (!p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.rightPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.leftPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							//System.err.println(z[pindex]+"\t"+lik.get(i-(relright/binsize))+"\t"+lik.get(j-(relleft/binsize)));
							tor[i][j] += z[pindex]*lik.get(i-(relright/binsize))*lik.get(j-(relleft/binsize));
						}
					}
				}
				if (!p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.rightPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.leftPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relright/binsize)-i)*lik.get(j-(relleft/binsize));
						}
					}
				}
				if (p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.rightPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.leftPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relright/binsize))*lik.get((relleft/binsize)-j);
						}
					}
				}
				if (p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.rightPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.rightPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.leftPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.leftPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relright/binsize)-i)*lik.get((relleft/binsize)-j);
						}
					}
				}
				goodset.add(p);
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			if (e1.contains(p.leftPointRegion(g)) && !goodset.contains(p)) {
				int pindex = index.get(p);
				int relright = p.rightPos-r2.getStart();
				int relleft = p.leftPos-r1.getStart();
				if (!p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.leftPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.rightPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relleft/binsize))*lik.get(j-(relright/binsize));
						}
					}
				}
				if (!p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.leftPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos)-r2.getStart())/binsize; j<(Math.min(r2.getEnd(), p.rightPos+lik.size()*binsize)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relleft/binsize)-i)*lik.get(j-(relright/binsize));
						}
					}
				}
				if (p.leftStrand && !p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos)-r1.getStart())/binsize; i<(Math.min(r1.getEnd(),p.leftPos+lik.size()*binsize)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.rightPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get(i-(relleft/binsize))*lik.get((relright/binsize)-j);
						}
					}
				}
				if (p.leftStrand && p.rightStrand) {
					for (int i=(Math.max(r1.getStart(), p.leftPos-lik.size()*binsize)-r1.getStart())/binsize+1; i<(Math.min(r1.getEnd(),p.leftPos)-r1.getStart())/binsize; i++) {
						for (int j=(Math.max(r2.getStart(), p.rightPos-lik.size()*binsize)-r2.getStart())/binsize+1; j<(Math.min(r2.getEnd(), p.rightPos)-r2.getStart())/binsize; j++) {
							tor[i][j] += z[pindex]*lik.get((relleft/binsize)-i)*lik.get((relright/binsize)-j);
						}
					}
				}
			}
		}
		return tor;
	}

	public double[] getForwardProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(i-pmkern);
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(i-pmkern);
				}
			}
		}
		for (int i=0; i<tor.length; i++) {
			tor[i] /= ((double)posTotal);
		}
		return tor;
	}

	public double fracLessThan(Region r) {
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), 0, true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		double count = (double)leftset.subSet(from, to).size();
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), 0, true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0f);
		count += (double) rightset.subSet(from, to).size();
		return count / (2d*((double)leftset.size()));
	}

	public double[] getUnnormForwardProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(otherkernel.size(), otherkernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(i-pmkern);
				}
			} else {
				for (int i=Math.max(r.getStart(),p.leftPos-otherkernel.size()+1); i<Math.min(r.getEnd(), p.leftPos); i++) {
					tor[i-r.getStart()] += otherkernel.get(p.leftPos-i);
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(i-pmkern);
				}
			} else {
				for (int i=Math.max(r.getStart(),p.rightPos-otherkernel.size()+1); i<Math.min(r.getEnd(), p.rightPos); i++) {
					tor[i-r.getStart()] += otherkernel.get(p.rightPos-i);
				}
			}
		}
		return tor;
	}

	/*
	 * returns an unnormalized probability
	 */
	public double[] getForwardWeightedProfile(Region r, double[] z) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.leftStrand) {
				double weight = z[index.get(p)];
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += weight*kernel.get(i-pmkern);
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (!p.rightStrand) {
				if (index.get(p)==null) System.err.println(p);
				double weight = z[index.get(p)];
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += weight*kernel.get(i-pmkern);
				}
			}
		}
		return tor;
	}

	/*
	 * prob of seeing r given r2
	 * normalized given r2
	 */
	public double[] getForwardWeightedConditionalProfile(Region r, double[] z, Region r2, double regionwidth) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(halfkern, otherkernel.size());
		Region e2 = r2.expand(halfkern, otherkernel.size());
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			double mult;

			if (!p.leftStrand) {
				int tmp = p.leftPos-r2.getStart()+halfkern;
				if (tmp<0||tmp>=kernel.size()) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = p.leftPos-r2.getStart();
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.rightPointRegion(g))) {

				if (!p.rightStrand) {
					int pmkern = p.rightPos-halfkern;
					int ppkern = p.rightPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(i-pmkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.rightPos-otherkernel.size()); i<Math.min(r.getEnd(), p.rightPos); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(p.rightPos-i);
					}
				}

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double mult;
			if (!p.rightStrand) {
				int tmp = p.rightPos-r2.getStart()+halfkern;
				if (tmp<0||tmp>=kernel.size()) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = p.rightPos-r2.getStart();
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.leftPointRegion(g))) {

				if (!p.leftStrand) {
					int pmkern = p.leftPos-halfkern;
					int ppkern = p.leftPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(i-pmkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.leftPos-otherkernel.size()); i<Math.min(r.getEnd(), p.leftPos); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(p.leftPos-i);
					}
				}

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		return tor;
	}

	/*
	 * prob of seeing r given r2
	 * normalized given r2
	 */
	public double[] getReverseForwardWeightedConditionalProfile(Region r, double[] z, Region r2, double regionwidth) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(halfkern, otherkernel.size());
		Region e2 = r2.expand(otherkernel.size()-1, halfkern);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			double mult;
			if (p.leftStrand) {
				int tmp = p.leftPos-r2.getStart()+halfkern;
				if (tmp>=kernel.size()||tmp<0) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = r2.getStart()-p.leftPos;
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.rightPointRegion(g))) {

				if (!p.rightStrand) {
					int pmkern = p.rightPos-halfkern;
					int ppkern = p.rightPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(i-pmkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.rightPos-otherkernel.size()); i<Math.min(r.getEnd(), p.rightPos); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(p.rightPos-i);
					}
				}

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double mult;
			if (p.rightStrand) {
				int tmp = p.rightPos-r2.getStart()+halfkern;
				if (tmp>=kernel.size()||tmp<0) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = r2.getStart()-p.rightPos;
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.leftPointRegion(g))) {

				if (!p.leftStrand) {
					int pmkern = p.leftPos-halfkern;
					int ppkern = p.leftPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(i-pmkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.leftPos-otherkernel.size()); i<Math.min(r.getEnd(), p.leftPos); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(p.leftPos-i);
					}
				}

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		return tor;
	}

	public double[] getReverseProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(kernel.size()-1-(i-pmkern));
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(kernel.size()-1-(i-pmkern));
				}
			}
		}
		for (int i=0; i<tor.length; i++) {
			tor[i] /= ((double)negTotal);
		}
		return tor;
	}

	public double[] getUnnormReverseProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(otherkernel.size(), otherkernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.leftStrand) {
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(p.leftPos-i+halfkern);
				}
			} else {
				for (int i=Math.max(r.getStart(),p.leftPos); i<Math.min(r.getEnd(), p.leftPos+otherkernel.size()); i++) {
					tor[i-r.getStart()] += otherkernel.get(i-p.leftPos);
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.rightStrand) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += kernel.get(kernel.size()-1-(i-pmkern));
				}
			} else {
				for (int i=Math.max(r.getStart(),p.rightPos); i<Math.min(r.getEnd(), p.rightPos+otherkernel.size()); i++) {
					tor[i-r.getStart()] += otherkernel.get(i-p.rightPos);
				}
			}
		}
		return tor;
	}

	public double[] getUnnormSymProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size()/2, kernel.size()/2);
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		//System.err.println(subset.size()+" in subset");
		for (PairedHit p : subset) {
			int pmkern = p.leftPos-halfkern;
			int ppkern = p.leftPos+halfkern;
			for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
				tor[i-r.getStart()] += kernel.get(p.leftPos-i+halfkern);
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		//System.err.println(subset.size()+" in subset");
		for (PairedHit p : subset) {
			int pmkern = p.rightPos-halfkern;
			int ppkern = p.rightPos+halfkern;
			for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
				tor[i-r.getStart()] += kernel.get(p.rightPos-i+halfkern);
			}
		}
		return tor;
	}
	
	public double[] getUnnormSymProfile(Region r, int window) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(window/2, window/2);
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		//System.err.println(subset.size()+" in subset");
		double inwin = 1d / ((double)window);
		for (PairedHit p : subset) {
			int pmkern = p.leftPos-(window/2);
			int ppkern = p.leftPos+(window/2);
			for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
				tor[i-r.getStart()] += inwin;
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		//System.err.println(subset.size()+" in subset");
		for (PairedHit p : subset) {
			int pmkern = p.rightPos-(window/2);
			int ppkern = p.rightPos+(window/2);
			for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
				tor[i-r.getStart()] += inwin;
			}
		}
		return tor;
	}

	public double[] getReverseWeightedProfile(Region r, double[] z) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.leftStrand) {
				double weight = z[index.get(p)];
				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += weight*kernel.get(kernel.size()-1-(i-pmkern));
				}
			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			if (p.rightStrand) {
				double weight = z[index.get(p)];
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					tor[i-r.getStart()] += weight*kernel.get(kernel.size()-1-(i-pmkern));
				}
			}
		}
		for (int i=0; i<tor.length; i++) {
			tor[i] /= ((double)negTotal);
		}
		return tor;
	}

	/*
	 * prob of seeing r given r2
	 * normalized given r2
	 */
	public double[] getReverseWeightedConditionalProfile(Region r, double[] z, Region r2, double regionwidth) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(otherkernel.size(), halfkern);
		Region e2 = r2.expand(halfkern, otherkernel.size());
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			double mult;
			if (!p.leftStrand) {
				int tmp = p.leftPos-r2.getStart()+halfkern;
				if (tmp<0||tmp>=kernel.size()) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = p.leftPos-r2.getStart();
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];
			norm += mult*weight;
			if (e.contains(p.rightPointRegion(g))) {
				if (p.rightStrand) {
					int pmkern = p.rightPos-halfkern;
					int ppkern = p.rightPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(p.rightPos-i+halfkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.rightPos); i<Math.min(r.getEnd(), p.rightPos+otherkernel.size()); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(i-p.rightPos);
					}
				}

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double mult;
			if (!p.rightStrand) {
				int tmp = p.rightPos-r2.getStart()+halfkern;
				if (tmp<0||tmp>=kernel.size()) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = p.rightPos-r2.getStart();
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.leftPointRegion(g))) {

				if (p.leftStrand) {
					int pmkern = p.leftPos-halfkern;
					int ppkern = p.leftPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(p.leftPos-i+halfkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.leftPos); i<Math.min(r.getEnd(), p.leftPos+otherkernel.size()); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(i-p.leftPos);
					}
				}

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		return tor;
	}

	/*
	 * prob of seeing r given r2
	 * normalized given r2
	 */
	public double[] getReverseReverseWeightedConditionalProfile(Region r, double[] z, Region r2, double regionwidth) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(otherkernel.size(), halfkern);
		Region e2 = r2.expand(otherkernel.size()-1, halfkern);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			double mult;
			if (p.leftStrand) {
				int tmp = p.leftPos-r2.getStart()+halfkern;
				if (tmp>=kernel.size()||tmp<0) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = r2.getStart()-p.leftPos;
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];
			norm += mult*weight;
			if (e.contains(p.rightPointRegion(g))) {
				if (p.rightStrand) {
					int pmkern = p.rightPos-halfkern;
					int ppkern = p.rightPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(p.rightPos-i+halfkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.rightPos); i<Math.min(r.getEnd(), p.rightPos+otherkernel.size()); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(i-p.rightPos);
					}
				}

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double mult;
			if (p.rightStrand) {
				int tmp = p.rightPos-r2.getStart()+halfkern;
				if (tmp>=kernel.size()||tmp<0) {
					mult = 0;
				} else {
					mult = kernel.get(tmp);
				}
			} else {
				int tmp = r2.getStart()-p.rightPos;
				if (tmp<0) {
					mult = 0;
				} else {
					mult = otherkernel.get(tmp);
				}
			}
			double weight = z[index.get(p)];

			norm += mult*weight;

			if (e.contains(p.leftPointRegion(g))) {

				if (p.leftStrand) {
					int pmkern = p.leftPos-halfkern;
					int ppkern = p.leftPos+halfkern;
					for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
						tor[i-r.getStart()] += mult*weight*kernel.get(p.leftPos-i+halfkern);
					}
				} else {
					for (int i=Math.max(r.getStart(),p.leftPos); i<Math.min(r.getEnd(), p.leftPos+otherkernel.size()); i++) {
						tor[i-r.getStart()] += mult*weight*otherkernel.get(i-p.leftPos);
					}
				}

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		return tor;
	}

	public double[] getSymWeightedConditionalProfile(Region r, double[] z, Region r2) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(halfkern, halfkern);
		Region e2 = r2.expand(halfkern, halfkern);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		for (PairedHit p : subset) {
			double leftdist = Math.abs(p.leftPos-r2.getStart());
			double weight = z[index.get(p)];
			norm += weight;
			if (e.contains(p.rightPointRegion(g))) {
				int pmkern = p.rightPos-halfkern;
				int ppkern = p.rightPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					double rightdist = Math.abs(i-p.rightPos);
					int threedist = (int)Math.floor(Math.sqrt(leftdist*leftdist+rightdist*rightdist));
					if (threedist<=halfkern) {
						tor[i-r.getStart()] += weight*kernel.get(threedist+halfkern);
					}
				}

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double rightdist = Math.abs(p.rightPos-r2.getStart());
			double weight = z[index.get(p)];

			norm += weight;

			if (e.contains(p.leftPointRegion(g))) {

				int pmkern = p.leftPos-halfkern;
				int ppkern = p.leftPos+halfkern;
				for (int i=Math.max(r.getStart(), pmkern); i<Math.min(r.getEnd(),ppkern); i++) {
					double leftdist = Math.abs(i-p.leftPos);
					int threedist = (int)Math.floor(Math.sqrt(leftdist*leftdist+rightdist*rightdist));
					if (threedist<=halfkern) {
						tor[i-r.getStart()] += weight*kernel.get(threedist+halfkern);
					}
				}

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		return tor;
	}

	public double getTmpNorm() {
		return tmpnorm;
	}
	
	public double getSymWeightedConditionalProfile(Region r, double[] z, Region r2, int window) {
		double[] tor = new double[r.getWidth()];
		int halfwin = window/2;
		Region e = r.expand(halfwin, halfwin);
		Region e2 = r2.expand(halfwin, halfwin);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		double numer = 0;
		double inwin = 1d / ((double)window);
		for (PairedHit p : subset) {
			double leftdist = Math.abs(p.leftPos-r2.getStart());
			double weight = z[index.get(p)];
			norm += weight;
			if (e.contains(p.rightPointRegion(g))) {
				numer += inwin*weight;

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double rightdist = Math.abs(p.rightPos-r2.getStart());
			double weight = z[index.get(p)];

			norm += weight;

			if (e.contains(p.leftPointRegion(g))) {

				numer += inwin*weight;

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		//return numer/norm;
		return numer;
	}
	
	public double getSymWeightedConditionalProfile(Region r, Region r2, int window) {
		double[] tor = new double[r.getWidth()];
		int halfwin = window/2;
		Region e = r.expand(halfwin, halfwin);
		Region e2 = r2.expand(halfwin, halfwin);
		PairedHit from = new PairedHit(g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		double norm = 0;
		double numer = 0;
		double inwin = 1d / ((double)window);
		for (PairedHit p : subset) {
			double leftdist = Math.abs(p.leftPos-r2.getStart());
			double weight = 1d;
			norm += weight;
			if (e.contains(p.rightPointRegion(g))) {
				numer += inwin*weight;

			}
		}
		from = new PairedHit(0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getStart(), true, (short)0,  0f);
		to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r2.getChrom()), e2.getEnd(), true, (short)0, 0f);
		subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			double rightdist = Math.abs(p.rightPos-r2.getStart());
			double weight = 1d;

			norm += weight;

			if (e.contains(p.leftPointRegion(g))) {

				numer += inwin*weight;

			}
		}
		/*
		for (int i=0; i<tor.length; i++) {
			tor[i] = (tor[i]+norm/regionwidth)/(2d*norm);
		}
		 */
		tmpnorm = norm;
		//return numer/norm;
		return numer;
	}

	/*
	public double[] getLeftProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(g.getChromID(r.getChrom()), e.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit to = new PairedHit(g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		SortedSet<PairedHit> subset = leftset.subSet(from, to);
		for (PairedHit p : subset) {
			int pmkern = p.leftPos-halfkern;
			int ppkern = p.leftPos+halfkern;
			for (int i=Math.max(e.getStart(), pmkern); i<Math.min(e.getEnd(),ppkern); i++) {
				tor[i-e.getStart()] += kernel.get(i-ppkern);
			}
		}
		return tor;
	}

	public double[] getRightProfile(Region r) {
		double[] tor = new double[r.getWidth()];
		Region e = r.expand(kernel.size(), kernel.size());
		PairedHit from = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getStart(), true, (short)0,  0f);
		PairedHit to = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), e.getEnd(), true, (short)0, 0f);
		SortedSet<PairedHit> subset = rightset.subSet(from, to);
		for (PairedHit p : subset) {
			int pmkern = p.rightPos-halfkern;
			int ppkern = p.rightPos+halfkern;
			for (int i=Math.max(e.getStart(), pmkern); i<Math.min(e.getEnd(),ppkern); i++) {
				tor[i-e.getStart()] += kernel.get(i-ppkern);
			}
		}
		return tor;
	}
	 */

	public int getMarginalCount(Region r) {
		PairedHit leftfrom = new PairedHit(g.getChromID(r.getChrom()), r.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(r.getChrom()), r.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(r.getChrom()), r.getEnd(), true, (short)0, 0f);
		Set<PairedHit> pointset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit hit : subset) {
			pointset.add(hit);
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit hit : subset) {
			pointset.add(hit);
		}
		return pointset.size();
	}

	public int getCount(Region anchor, Region region) {
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		Set<PairedHit> pointset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
				//pointset.add(new Point(g,region.getChrom(),p.rightPos));
				pointset.add(p);
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
				//pointset.add(new Point(g,region.getChrom(),p.leftPos));
				pointset.add(p);
			}
		}
		return pointset.size();
	}
	
	public Set<PairedHit> getHitSet(Region anchor, Region region) {
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		Set<PairedHit> pointset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
				//pointset.add(new Point(g,region.getChrom(),p.rightPos));
				pointset.add(p);
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
				//pointset.add(new Point(g,region.getChrom(),p.leftPos));
				pointset.add(p);
			}
		}
		return pointset;
	}
	
	public double computeRatio(Region anchor, Region region) {
		Set<PairedHit> tmpset = getHitSet(anchor,region);
		double tor = 0;
		for (PairedHit hit : tmpset) {
			tor += fullMap.get(hit) - restrictedMap.get(hit);
		}
		return tor;
	}
	
	public double computeRatio(Region anchor, Region region, double modifier) {
		Set<PairedHit> tmpset = getHitSet(anchor,region);
		double tor = 0;
		for (PairedHit hit : tmpset) {
			tor += fullMap.get(hit) - restrictedMap.get(hit) - modifier;
		}
		return tor;
	}
	
	public double fullValue(PairedHit hit) {
		return fullMap.get(hit);
	}
	
	public double restrictedValue(PairedHit hit) {
		return restrictedMap.get(hit);
	}

	public double getCount(Region anchor, Region region, double[] z) {
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		Set<PairedHit> pointset = new HashSet<PairedHit>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
				//pointset.add(new Point(g,region.getChrom(),p.rightPos));
				pointset.add(p);
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
				//pointset.add(new Point(g,region.getChrom(),p.leftPos));
				pointset.add(p);
			}
		}
		double sum = 0;
		for (PairedHit p : pointset) {
			sum += z[index.get(p)];
		}
		return sum;
	}

	public int getCount(Region anchor, int anchorcenter, Region region) {
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		SortedSet<Point> pointset = new TreeSet<Point>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if ((p.leftStrand && p.leftPos>=anchorcenter) || (!p.leftStrand && p.leftPos<=anchorcenter)) {
				if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
					pointset.add(new Point(g,region.getChrom(),p.rightPos));
				}
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if ((p.rightStrand && p.rightPos>=anchorcenter) || (!p.rightStrand && p.rightPos<=anchorcenter)) {
				if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
					pointset.add(new Point(g,region.getChrom(),p.leftPos));
				}
			}
		}
		return pointset.size();
	}

	public SortedSet<Point> getStrandedPointSet(Region anchor, int anchorcenter, Region region) {
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		SortedSet<Point> pointset = new TreeSet<Point>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if ((p.leftStrand && p.leftPos>=anchorcenter) || (!p.leftStrand && p.leftPos<=anchorcenter)) {
				if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
					pointset.add(new Point(g,region.getChrom(),p.rightPos));
				}
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if ((p.rightStrand && p.rightPos>=anchorcenter) || (!p.rightStrand && p.rightPos<=anchorcenter)) {
				if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
					pointset.add(new Point(g,region.getChrom(),p.leftPos));
				}
			}
		}
		return pointset;
	}

	public Map<Region,Double> test(Region anchor, Region region) {
		int anchorchrom = anchor.getGenome().getChromID(anchor.getChrom());
		int regionchrom = region.getGenome().getChromID(region.getChrom());
		PairedHit leftfrom = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit leftto = new PairedHit(g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0, 0, true, (short)0, 0f);
		PairedHit rightfrom = new PairedHit(0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getStart(), true, (short)0,  0f);
		PairedHit rightto = new PairedHit( 0, 0, true, (short)0, g.getChromID(anchor.getChrom()), anchor.getEnd(), true, (short)0, 0f);
		SortedSet<Point> pointset = new TreeSet<Point>();
		SortedSet<PairedHit> subset = leftset.subSet(leftfrom, leftto);
		for (PairedHit p : subset) {
			if (p.rightChrom==regionchrom && region.getStart()<=p.rightPos && region.getEnd()>=p.rightPos) {
				pointset.add(new Point(g,region.getChrom(),p.rightPos));
			}
		}
		subset = rightset.subSet(rightfrom, rightto);
		for (PairedHit p : subset) {
			if (p.leftChrom==regionchrom && region.getStart()<=p.leftPos && region.getEnd()>=p.leftPos) {
				pointset.add(new Point(g,region.getChrom(),p.leftPos));
			}
		}
		Map<Region,Double> pvals = new HashMap<Region,Double>();

		if (pointset.size()>0) {
			Iterator<Point> iter = pointset.iterator();
			Point first = iter.next();
			Region probe = new Region(g, region.getChrom(), first.getLocation(), first.getLocation()+1);
			LinkedList<Point> probeset = new LinkedList<Point>();
			probeset.addFirst(first);
			while (iter.hasNext()) {
				while (iter.hasNext() && probe.getWidth()<minsize) {
					Point next = iter.next();
					probeset.addLast(next);
					Region nextr = new Region(g, region.getChrom(), next.getLocation(), next.getLocation()+1);
					probe = probe.combine(nextr);
				}

				while (probe.getWidth()>maxsize) {
					probeset.removeFirst();
					probe = new Region(g,probe.getChrom(),probeset.getFirst().getLocation(),probeset.getLast().getLocation());
				}

				if (probe.getWidth()>=minsize) {
					double pval = pval(anchor, probe, probeset.size());
					pvals.put(probe,pval);
					probeset.removeFirst();
					probe = new Region(g,probe.getChrom(),probeset.getFirst().getLocation(),probeset.getLast().getLocation());
				}
			}
		}
		return pvals;
	}

	public Map<Region,Double> test(Region anchor, int anchorcenter, Region region) {
		SortedSet<Point> pointset = getStrandedPointSet(anchor, anchorcenter, region);
		Map<Region,Double> pvals = new HashMap<Region,Double>();

		if (pointset.size()>0) {
			Iterator<Point> iter = pointset.iterator();
			Point first = iter.next();
			Region probe = new Region(g, region.getChrom(), first.getLocation(), first.getLocation()+1);
			LinkedList<Point> probeset = new LinkedList<Point>();
			probeset.addFirst(first);
			while (iter.hasNext()) {
				while (iter.hasNext() && probe.getWidth()<minsize) {
					Point next = iter.next();
					probeset.addLast(next);
					Region nextr = new Region(g, region.getChrom(), next.getLocation(), next.getLocation()+1);
					probe = probe.combine(nextr);
				}

				while (probe.getWidth()>maxsize) {
					probeset.removeFirst();
					probe = new Region(g,probe.getChrom(),probeset.getFirst().getLocation(),probeset.getLast().getLocation());
				}

				if (probe.getWidth()>=minsize) {
					double pval = pval(anchor, anchorcenter, probe, probeset.size());
					pvals.put(probe,pval);
					probeset.removeFirst();
					probe = new Region(g,probe.getChrom(),probeset.getFirst().getLocation(),probeset.getLast().getLocation());
				}
			}
		}
		return pvals;
	}

	private double sum(double[] arr) {
		double sum = 0;
		for (int i=0; i<arr.length; i++) {
			sum += arr[i];
		}
		return sum;
	}

	public double pval(Region anchor, Region probe, int count) {
		double[] leftanchorprofile = getForwardProfile(anchor);
		double[] rightanchorprofile = getReverseProfile(anchor);
		double prob1 = (sum(leftanchorprofile) + sum(rightanchorprofile))/2.0d;
		double[] leftprobeprofile = getForwardProfile(probe);
		double[] rightprobeprofile = getReverseProfile(probe);
		double prob2 = (sum(leftprobeprofile) + sum(rightprobeprofile))/2.0d;

		poisson.setMean(chimericReads*prob1*prob2);
		return 1.0d - poisson.cdf(count) + poisson.pdf(count);
	}

	public double pval(Region anchor, int anchorcenter, Region probe, int count) {
		double[] leftanchorprofile = getForwardProfile(new Region(g,anchor.getChrom(),anchorcenter,anchor.getEnd()));
		double[] rightanchorprofile = getReverseProfile(new Region(g,anchor.getChrom(),anchor.getStart(),anchorcenter));
		double prob1 = (sum(leftanchorprofile) + sum(rightanchorprofile))/2.0d;
		double[] leftprobeprofile = getForwardProfile(probe);
		double[] rightprobeprofile = getReverseProfile(probe);
		double prob2 = (sum(leftprobeprofile) + sum(rightprobeprofile))/2.0d;

		poisson.setMean(chimericReads*prob1*prob2);
		return 1.0d - poisson.cdf(count) + poisson.pdf(count);
	}

}
