package edu.mit.csail.cgs.metagenes;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.StrandedPoint;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.verbs.GenomeExpander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.RefGeneGenerator;

//Only contains some loaders right now

public class MetaUtils {

	private Genome genome;
	
	public MetaUtils(Genome g){
		genome = g;
	}
	
	public Iterator<Point> loadTSSs() { 
		// TODO: By default, loads the refseq annotations for this genome.  But this won't work 
		// for all genomes, therefore, this is a bit of a hack. We should update this to be 
		// more general...
		System.out.println("Loading gene TSSs");
		RefGeneGenerator gen = new RefGeneGenerator(genome, "refGene");
		Mapper<Gene,Point> tssMapper = new Mapper<Gene,Point>() { 
			public Point execute(Gene g) { 
				if(g.getStrand() == '+') { 
					return new StrandedPoint(g.getGenome(), g.getChrom(), g.getStart(), '+');
				} else { 
					return new StrandedPoint(g.getGenome(), g.getChrom(), g.getEnd(), '-');
				}
			}
		};
		GenomeExpander<Gene> gexp = new GenomeExpander<Gene>(gen);
		Iterator<Gene> genes = gexp.execute(genome);
		Iterator<Point> points = new MapperIterator<Gene,Point>(tssMapper, genes);
		return points;
	}
	
	public Vector<Point> loadPoints(File f) throws IOException {
		System.out.println("Loading points");
		Vector<Point> pts = new Vector<Point>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null) {
			if (line.trim().equals(""))
				continue;
			String [] curr = line.split("\\s+");
			String coord = curr[0];
			if(curr.length>=3 && curr[2].contains(":")){coord = curr[2];}
			
			if(coord.contains(":")) {
				String [] currB = coord.split(":");
				String chrom = currB[0];
				char strand = '?';
				if(currB.length==3)
					strand = currB[2].charAt(0);
				Point pt = null;
				int location=-1;
				if(currB[1].contains("-")){
					String [] currC = currB[1].split("-");
					location = (new Integer(currC[0])+new Integer(currC[1]))/2; 
				}else{
					location = new Integer(currB[1]);
				}
				if(strand!='?')
					pt = new StrandedPoint(genome, chrom, location, strand);
				else
					pt = new Point(genome, chrom, location);
			
				pts.add(pt);
			} else { 
				System.err.println(String.format("Couldn't find point in line \"%s\"", line));
			}
		}
		br.close();
		System.err.println(pts.size()+" points loaded");
		return pts;
	}
}
