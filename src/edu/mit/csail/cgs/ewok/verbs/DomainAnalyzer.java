package edu.mit.csail.cgs.ewok.verbs;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingScan;
import edu.mit.csail.cgs.datasets.binding.BindingScanLoader;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.chipchip.SQLGeneric;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.GeneDomainData;
import edu.mit.csail.cgs.ewok.nouns.GeneDomainTimeSeries;
import edu.mit.csail.cgs.ewok.nouns.SimpleDomain;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.NotFoundException;

public class DomainAnalyzer implements Closeable {
	
	public static String[][] mark_stages = { 
		{ "H3K27me3", "ES Stage", "es" },
		{ "H3K27me3", "ES+2d Stage, before RA", "es2" },
		{ "H3K27me3", "Olig2 Stage", "olig2" },
		{ "H3K27me3", "Hb9 Stage", "hb9" },

		{ "H3K4me3", "ES Stage", "es" },
		{ "H3K4me3", "ES+2d Stage, before RA", "es2" },
		{ "H3K4me3", "Olig2 Stage", "olig2" },
		{ "H3K4me3", "Hb9 Stage", "hb9" },
	};
	
	public static void main(String[] args) { 
		try {

			Genome mm8 = Organism.findGenome("mm8");
			Collection<Gene> genes = null;
			LinkedList<GeneDomainTimeSeries> tss = 
				new LinkedList<GeneDomainTimeSeries>();

			for(int i = 0; i < mark_stages.length; i++) {
				String mark = mark_stages[i][0], stage = mark_stages[i][1];
				String stageKey = mark_stages[i][2];
				File output = new File(mark + "_" + stageKey + "_domains.txt");

				System.out.println("Processing " + mark + "," + stage);
				PrintStream ps = new PrintStream(new FileOutputStream(output));
				
				DomainAnalyzer da = new DomainAnalyzer(mm8, mark, stage, stageKey);
				if(genes == null) { 
					genes = da.getWellTiledGenes();
					for(Gene g : genes) { 
						tss.addLast(new GeneDomainTimeSeries(g, 8));
					}
				}
				
				//da.doAnalysis();
				
				Collection<SimpleDomain> doms = da.findAllDomains(ps);
				ps.println("# domains: " + doms.size());
				System.out.println("\t# Domains: " + doms.size());
				
				ps.close();
				
				da.enterAsBindingScan("SimpleDomain", doms);
				
				for(SimpleDomain dom : doms) { 
					for(GeneDomainTimeSeries gdts : tss) { 
						gdts.addDomain(dom, i);
					}
				}
			}
			
			PrintStream ps = new PrintStream(new FileOutputStream("time_series.txt"));
			for(GeneDomainTimeSeries gd : tss) { 
				ps.println(gd.toString());
			}
			ps.close();
			
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}
	
	public static int window = 30000;

	private Genome genome;
	private String mark, stage, stageKey;
	private ChipChipLocator loc;
	private ChipChipData data;
	private Expander<Region,SimpleDomain> domainCaller;
	private RefGeneGenerator<Region> geneGen;
	
	public DomainAnalyzer(Genome g, String m, String s, String sk) { 
		genome = g;
		mark = m; stage = s; stageKey = sk;
		String name = "Mm " + mark + ":HBG3:" + stage + 
			" vs H3:HBG3:" + stage;
		String version = "median linefit";

		loc = new ChipChipLocator(genome, name, version);
		data = loc.createObject();
		geneGen = new RefGeneGenerator<Region>(genome, "refGene");

		SimpleDomainFinder sdf = new SimpleDomainFinder(data);
		if(mark.equals("H3K27me3") && stage.startsWith("ES")) { 
			sdf.setMeanThreshold(2.0);
		}

		//domainCaller = sdf;
		domainCaller = new HMMDomainGenerator(loc);
	}
	
	public boolean isClosed() { 
		return data == null;
	}
	
	public void close() { 
		((SQLData)data).close();
		data = null;
	}
	
	public void enterAsBindingScan(String type, Collection<? extends Region> regions) 
		throws IOException, SQLException {
		
		String version = loc.name;

		if(regions == null) { 
			File f = new File(mark + "_" + stageKey + "_domains.txt");
			regions = loadRegions(f);
		}

		BindingScanLoader loader = new BindingScanLoader();
		Connection cxn = loader.getConnection();
		cxn.setAutoCommit(false);
		
		BindingScan scan = new BindingScan(genome, type, version);
			
		Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
		LinkedList<Region> scannedRegions = new LinkedList<Region>();
		while(chroms.hasNext()) { 
			scannedRegions.addLast(chroms.next());
		}
		
		Map<String,String> params = new HashMap<String,String>();
		params.put("author", "DomainAnalyzer");

		loader.insertScan(scan);
		loader.insertNewRegions(scan, scannedRegions);
		loader.insertNewParams(scan, params);
		//loader.insertExpt(scan, loc);
		
		int exptType = BindingScan.getLocatorType(loc);
		int[] exptIDs = BindingScan.getChipChipIDs(loc, cxn);
		loader.insertNewExpts(scan, exptIDs, exptType);
		
		for(Region r : regions) { 
			SimpleDomain sd = new SimpleDomain(r);
			BindingEvent evt = sd.getBindingEvent();
			loader.insertEvent(scan, evt);
		}

		System.out.println("Inserted " + scan.toString() + " (" + regions.size() + " regions)");

		cxn.commit();
		cxn.setAutoCommit(true);
		loader.close();
	}
	
	private LinkedList<Region> loadRegions(File f) throws IOException { 
		LinkedList<Region> regions = new LinkedList<Region>();
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = null;
		while((line = br.readLine()) != null) { 
			line = line.trim();
			if(line.length() > 0 && !line.startsWith("#")) { 
				String[] array = line.split("\\s+");
				String chrom = array[0];
				int start = Integer.parseInt(array[1]);
				int end = Integer.parseInt(array[2]);
				Region r = new Region(genome, chrom, start, end);
				regions.addLast(r);
			}
		}
		
		br.close();
		
		return regions;
	}
	
	public Collection<SimpleDomain> findAllDomains(PrintStream ps) { 
		LinkedList<SimpleDomain> doms = new LinkedList<SimpleDomain>();
		
		//ChromRegionIterator itr = new ChromRegionIterator(genome);
		//Iterator<Region> regions = new MapperIterator<NamedRegion,Region>(
				//new CastingMapper<NamedRegion,Region>(), itr);
		
		Iterator<Region> regions = new WellTiledRegionGenerator();

		Iterator<SimpleDomain> domains = 
			new ExpanderIterator<Region,SimpleDomain>(domainCaller, regions);
		
		int c = 0;
		while(domains.hasNext()) { 
			SimpleDomain dom = domains.next();
			doms.addLast(dom);
			c += 1;
			
			Iterator<Gene> genes = geneGen.execute(dom);
			int g = 0;
			while(genes.hasNext()) { genes.next(); g += 1; }
			
			ps.println(dom.getChrom() + "\t" + 
					dom.getStart() + "\t" + dom.getEnd() + "\t" + 
					dom.getWidth() + "\t" + g);
		}

		return doms;
	}
	
	public void doAnalysis() { 
		TiledGeneFilter tiledFilter = 
			new TiledGeneFilter(data, 10, window);
		
		DomainAnalysis<SimpleDomain> analysis = 
			new DomainAnalysis<SimpleDomain>(domainCaller, window);
		
		//Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
		//Iterator<Region> regions = new MapperIterator<NamedRegion,Region>( new CastingMapper<NamedRegion,Region>(), chroms);
		Iterator<Region> regions = new WellTiledRegionGenerator();

		Iterator<Gene> genes = 
			new ExpanderIterator<Region,Gene>(
					new RefGeneGenerator(genome, "refGene"), regions);
		
		genes = new FilterIterator<Gene,Gene>(tiledFilter, genes);
		Iterator<GeneDomainData> ditr = 
			new MapperIterator<Gene,GeneDomainData>(analysis, genes);
		
		while(ditr.hasNext()) { 
			GeneDomainData data = ditr.next();
			if(data.getNumDomains() > 0) { 
				data.printData();
			}
		}
	}
	
	public Collection<Gene> getWellTiledGenes() { 
		Iterator<Region> tiled = new WellTiledRegionGenerator();
		Iterator<Gene> genes = new ExpanderIterator<Region,Gene>(geneGen, tiled);
		Set<Gene> totalGenes = new TreeSet<Gene>();
		while(genes.hasNext()) { totalGenes.add(genes.next()); }
		return totalGenes;
	}
}

class WellTiledRegionGenerator implements Iterator<Region> { 

	private Genome genome;
	private File file;

	private LinkedList<Region> regions;

	public WellTiledRegionGenerator() { 
		try { 
			genome = Organism.findGenome("mm8");
		} catch(NotFoundException nfe) { 
			throw new IllegalStateException();
		}
		file = new File("well_tiled_regions.txt");
		regions = new LinkedList<Region>();
		parse();
	}

	public WellTiledRegionGenerator(File f) { 
		try { 
			genome = Organism.findGenome("mm8");
		} catch(NotFoundException nfe) { 
			throw new IllegalStateException();
		}
		file = f;
		regions = new LinkedList<Region>();
		parse();
	}

	public void parse() { 
		try { 
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = null;
			while((line = br.readLine()) != null) { 
				line = line.trim();
				String [] pieces = line.split("\t");
				if(line.length() > 0) { 
				 	String chrom = pieces[0];
					int start = Integer.parseInt(pieces[1]);
					int end = Integer.parseInt(pieces[2]);
					Region r = new Region(genome, chrom, start, end);
					regions.addLast(r);
				}
			}

			br.close();
		} catch(IOException ie) { 
			ie.printStackTrace(System.err);
		}
	}

	public boolean hasNext() { return !regions.isEmpty(); }
	public Region next() { return regions.removeFirst(); }
	public void remove() { throw new UnsupportedOperationException(); }
}

class TiledGeneFilter implements Filter<Gene,Gene> {

	private ChipChipData data;
	private int numProbes;
	private int window;
	
	public TiledGeneFilter(ChipChipData ccd, int np, int w) { 
		data = ccd;
		numProbes = np;
		window = w;
	}
	
	public Gene execute(Gene a) {
		try {
			data.window(a.getChrom(), a.getStart()-window, a.getEnd()+window);
			return data.getCount() >= numProbes ? a : null;
		} catch (NotFoundException e) {
			e.printStackTrace();
			return null;
		}
	}
}
