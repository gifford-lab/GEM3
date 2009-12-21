package edu.mit.csail.cgs.tools.chipchip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.utils.NotFoundException;

public class GeneFileIterator implements Iterator<Gene> {

	Iterator<Gene> iter;
	
	public GeneFileIterator(File geneFile) {
		List<Gene> geneList = new LinkedList<Gene>();
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(geneFile));
			String line;
			String[] split;
			Map<String,Genome> genMap = new HashMap<String,Genome>();
			while((line = br.readLine()) != null) {
				split = line.split("\t");
				String[] genSplit = split[0].split(",");
				Genome tmpGen;
				if (genMap.containsKey(genSplit[1])) {
					tmpGen = genMap.get(genSplit[1]);
				} else {
					tmpGen = Organism.findGenome(genSplit[1]);
					genMap.put(genSplit[1], tmpGen);
				}
				geneList.add(new Gene(tmpGen, split[1], Integer.parseInt(split[2]),
						Integer.parseInt(split[3]), split[4], split[5], split[6].charAt(0), split[7],
						Integer.parseInt(split[8])));
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
		iter = geneList.iterator();
	}
	
	public boolean hasNext() {
		return iter.hasNext();
	}

	public Gene next() {
		return iter.next();
	}

	public void remove() {
		iter.remove();
	}

}
