package edu.mit.csail.cgs.tools.hypotheses.utils;

import java.sql.SQLException;
import java.util.*;

import edu.mit.csail.cgs.datasets.general.Factor;
import edu.mit.csail.cgs.datasets.general.MetadataLoader;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.nouns.HarbisonRegCodeProbe;
import edu.mit.csail.cgs.ewok.verbs.CastingMapper;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.HarbisonRegCodeProbeGenerator;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;

public class HarbisonFactorIterator implements Iterator<Factor> {

	private Iterator<Factor> factors;
	
	public HarbisonFactorIterator(Genome genome) throws SQLException { 
		LinkedList<Factor> fs = new LinkedList<Factor>();
		MetadataLoader loader = new MetadataLoader();

		Iterator<NamedRegion> chroms = new ChromRegionIterator(genome);
		Iterator<Region> chromRegions = 
			new MapperIterator<NamedRegion,Region>(new CastingMapper<NamedRegion,Region>(), chroms);
		
		
		HarbisonRegCodeProbeGenerator gen = new HarbisonRegCodeProbeGenerator(genome);
		Iterator<HarbisonRegCodeProbe> probes = 
			new ExpanderIterator<Region,HarbisonRegCodeProbe>(gen, chromRegions);
		
		Set<String> nfs = new TreeSet<String>();

		while(probes.hasNext()) { 
			HarbisonRegCodeProbe probe = probes.next();
			for(String f : probe.getFactors()) { 
				if(!nfs.contains(f)) { 
					nfs.add(f);
				}
			}
		}
		
		for(String f : nfs) { 
			fs.addLast(loader.getFactor(f));
		}

		loader.close();
		factors = fs.iterator();
	}
	
	
	public boolean hasNext() {
		return factors.hasNext();
	}

	public Factor next() {
		return factors.next();
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}
