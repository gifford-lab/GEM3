package edu.mit.csail.cgs.tools.chipchip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.CoveredFilter;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.ExpanderIterator;
import edu.mit.csail.cgs.ewok.verbs.FilterIterator;
import edu.mit.csail.cgs.ewok.verbs.GeneToPromoter;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;
import edu.mit.csail.cgs.ewok.verbs.UniqueishGeneFilter;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;

public class CoveredGenes {

	private static final String[] geneTables = {"refGene", "ensGene", "mgcGenes"};
	
	/**
	 * @param args
	 * @throws NotFoundException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NotFoundException, IOException {
		Genome g = Args.parseGenome(args).cdr();
		int startBuffer = Args.parseInteger(args, "startbuf", 5);
		int endBuffer = Args.parseInteger(args, "endbuf", 5);
		int startProm = Args.parseInteger(args, "startprom", 500);
		int endProm = Args.parseInteger(args, "endprom", 500);
		String regionFile = Args.parseString(args, "regionfile", "");
		String outFile = Args.parseString(args, "outfile", "");
		String geneFile = Args.parseString(args, "genefile", "");
		System.err.println("read args");
		List<Gene> geneList = new LinkedList<Gene>();
		/*
		//PrintStream out = new PrintStream(outFile); //
		for (int i=0; i<geneTables.length; i++) {
			Expander<NamedRegion,Gene> gen = new RefGeneGenerator<NamedRegion>(g, geneTables[i]);
			Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
			Iterator<Gene> geneIter = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
			while (geneIter.hasNext()) {
				Gene tmp = geneIter.next();
				geneList.add(tmp);
				//out.println(tmp.getGenome()+"\t"+tmp.getChrom()+"\t"+tmp.getStart()+"\t"+tmp.getEnd()+"\t"+
				//		tmp.getName()+"\t"+tmp.getID()+"\t"+tmp.getStrand()+"\t"+tmp.getSource()+"\t"+
				//		tmp.getDBID()); //
			}
		}
		//out.flush(); //
		//out.close(); //
		 */
		BufferedReader br = new BufferedReader(new FileReader(geneFile));
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
		
		Iterator<Gene> geneIter = geneList.iterator();
		geneIter = new FilterIterator<Gene,Gene>(new UniqueishGeneFilter<Gene>(startBuffer,endBuffer),geneIter);
		Iterator<NamedStrandedRegion> promIter = new MapperIterator<Gene,NamedStrandedRegion>(new GeneToPromoter(startProm,endProm),geneIter);
		promIter = new FilterIterator<NamedStrandedRegion,NamedStrandedRegion>(new CoveredFilter<NamedStrandedRegion>(new File(regionFile)),promIter);
		
		PrintStream out = new PrintStream(outFile);
		while (promIter.hasNext()) {
			NamedStrandedRegion tmpProm = promIter.next();
			if (tmpProm!=null) {
				String name = tmpProm.getName().split(" ")[0];
				out.println(name+"\t"+tmpProm.getChrom()+"\t"+((tmpProm.getStart()+
						tmpProm.getEnd())/2)+"\t"+tmpProm.getStrand());
			}
		}
		out.flush();
		out.close();
		
	}

}
