package edu.mit.csail.cgs.ewok.verbs.sequence;

import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.utils.Coverage;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class RegionCoverage extends Coverage {

	private Genome genome;
	
	public RegionCoverage(Genome g) { 
		super();
		genome = g;
	}
	
	public RegionCoverage(RegionCoverageModel m) { 
		this(null, m);
	}
	
	public RegionCoverage(Genome g, RegionCoverageModel m) { 
		super(m);
		genome = g == null ? m.loadGenome() : g;
	}
	
	public RegionCoverage(Genome g, Coverage c) { 
		super(c);
		genome = g;
	}
	
	public RegionCoverage(Genome g, Collection<? extends Region> rs) { 
		this(g, rs, null);
	}
	
	public RegionCoverage(Genome g, Collection<? extends Region> rs, PrintStream ps) { 
		this(g);
		int i = 0;
		for(Region r : rs) { 
			addRegion(r);
			i += 1;
			if(ps != null) { 
				if(i % 100 == 0) { 
					System.out.print("."); System.out.flush();
				}
				if(i % 1000 == 0) { 
					System.out.println(String.format("(%dk)", i/1000));
				}
			}
		}
		if(ps != null) { ps.println(); }
	}
	
	public RegionCoverage(Genome g, Iterator<? extends Region> rs) { 
		this(g, rs, null);
	}
	
	public RegionCoverage(Genome g, Iterator<? extends Region> rs, PrintStream ps) { 
		this(g);
		int i = 0;
		while(rs.hasNext()) { 
			addRegion(rs.next());
			i += 1;
			if(ps != null) { 
				if(i % 100 == 0) { 
					System.out.print("."); System.out.flush();
				}
				if(i % 1000 == 0) { 
					System.out.println(String.format("(%dk)", i/1000));
				}
			}
		}
		
		if(ps != null) { ps.println(); }
	}
	
	public Iterator<Region> regions() { 
		Expander<String,Region> regExp = new Expander<String,Region>() { 
			public Iterator<Region> execute(String key) { 
				Iterator<Integer[]> units = units(key);
				return new MapperIterator<Integer[],Region>(new UnitMapper(key), units); 
			}
		};
		return new ExpanderIterator<String,Region>(regExp, keys());
	}
	
	private class UnitMapper implements Mapper<Integer[],Region> {
		private String key; 
		public UnitMapper(String k) { key = k; }
		public Region execute(Integer[] a) {
			return new Region(genome, key, a[0], a[1]);
		} 
	}
	
	public Region rightNearest(Point p) { 
		Integer[] bounds = super.rightNearest(p.getChrom(), p.getLocation());
		if(bounds != null) { 
			return new Region(p.getGenome(), p.getChrom(), bounds[0], bounds[1]);
		} else { 
			return null;
		}
	}
	
	public Region leftNearest(Point p) { 
		Integer[] bounds = super.leftNearest(p.getChrom(), p.getLocation());
		if(bounds != null) { 
			return new Region(p.getGenome(), p.getChrom(), bounds[0], bounds[1]);
		} else { 
			return null;
		}
	}
	
	public RegionCoverage subtract(Coverage c) { 
		return new RegionCoverage(genome, super.subtract(c));
	}
	
	public RegionCoverage intersection(Coverage c) { 
		return new RegionCoverage(genome, super.intersection(c));
	}
	
	public RegionCoverage union(Coverage c) { 
		return new RegionCoverage(genome, super.union(c));
	}
	
	public void addRegion(Region r) {
		if(!r.getGenome().equals(genome)) {
			throw new IllegalArgumentException();
		}
		addInterval(r.getChrom(), r.getStart(), r.getEnd());
	}
	
	public int coverage(Region r) { 
		return super.coverage(r.getChrom(), r.getStart(), r.getEnd());
	}
	
	public Collection<Region> coveredRegions(Region r) { 
		Collection<Integer[]> cov = covered(r.getChrom(), r.getStart(), r.getEnd());
		ArrayList<Region> regs = new ArrayList<Region>();
		for(Integer[] c : cov) { 
			regs.add(new Region(genome, r.getChrom(), c[0], c[1]));
		}
		return regs;
	}

	public boolean hasCoverage(Region r) {
		return hasCoverage(r.getChrom(), r.getStart(), r.getEnd());
	}
	
	public RegionCoverageModel asModel() { 
		return new RegionCoverageModel(this, super.asModel());
	}
	
	public static class RegionCoverageModel extends CoverageModel {
		
		public String species, version;
		
		public RegionCoverageModel() {}
		
		public RegionCoverageModel(RegionCoverage c, CoverageModel m) {
			keys = m.keys.clone();
			coverages = m.coverages.clone();
			species = c.genome.getSpecies();
			version = c.genome.getVersion();
		}
		
		public Genome loadGenome() { 
			try {
				Organism org = Organism.getOrganism(species);
				return org.getGenome(version);
			} catch (NotFoundException e) {
				e.printStackTrace();
				return null;
			}
		}
	}

	public boolean isContained(StrandedRegion r) {
		return isContained(r.getChrom(), r.getStart(), r.getEnd());
	}

	public Region findNearestRight(Point p) {
		Integer[] bounds = super.findNearestRight(p.getChrom(), p.getLocation());
		if(bounds == null) { return null; }
		return new Region(p.getGenome(), p.getChrom(), bounds[0], bounds[1]);
	}

	public Region findNearestLeft(Point p) {
		Integer[] bounds = super.findNearestLeft(p.getChrom(), p.getLocation());
		if(bounds == null) { return null; }
		return new Region(p.getGenome(), p.getChrom(), bounds[0], bounds[1]);
	}

	public Collection<Pair<Region, Region>> 
		findPairs(int maxDist, RegionCoverage coverage, boolean right) {

		ArrayList<Pair<Region,Region>> lst = new ArrayList<Pair<Region,Region>>();
		Iterator<String> chroms = keys();

		while(chroms.hasNext()) { 
			String chrom = chroms.next();
			Collection<Pair<Integer[],Integer[]>> 
				pairs = right ? super.findRightPairs(chrom, maxDist, coverage) 
						: super.findLeftPairs(chrom, maxDist, coverage); 
			
			for(Pair<Integer[],Integer[]> p : pairs) { 
				Integer[] p1 = p.getFirst(), p2 = p.getLast();
				System.out.println(String.format("%d,%d vs. %d,%d", p1[0], p1[1], p2[0], p2[1]));
				
				Region r1 = new Region(genome, chrom, p1[0], p1[1]);
				Region r2 = new Region(genome, chrom, p2[0], p2[1]);
				lst.add(new Pair<Region,Region>(r1, r2));
			}
		}
		
		return lst;
	}
}
