package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Write fasta sequence file
 * Input: a list of coordinates and a window size, or
 * 		  a list of regions
 * @author yuchun
 *
 */
public class FastaWriter{  
	// --species "Homo sapiens;hg19" --coords /your/path [ --win 100 --no_cache --use_db_genome --out name --genome genome_path ]
  public static void main(String[] args){
    ArgParser ap = new ArgParser(args);
    Set<String> flags = Args.parseFlags(args);
    Genome genome=null;
    
    try {
      Pair<Organism, Genome> pair = Args.parseGenome(args);
      if(pair==null){
        //Make fake genome... chr lengths provided???
        if(ap.hasKey("g")){
        	String gname = ap.getKeyValue("g").replaceFirst(".info", "");
          genome = new Genome(gname, new File(ap.getKeyValue("g")), true);
        }else{
              System.err.println("No genome provided; provide a Gifford lab DB genome name or a file containing chromosome name/length pairs.");;System.exit(1);
        }
      }else{
        genome = pair.cdr();
      }  
    } catch (NotFoundException e) {
      e.printStackTrace();
    }
    
//    boolean wantHMS = flags.contains("hms");
//    boolean wantChIPMunk = flags.contains("chipmunk");
    
    int window = Args.parseInteger(args, "win", 100);
    
    SequenceGenerator<Region> seqgen = new SequenceGenerator<Region>();
	seqgen.useCache(!flags.contains("no_cache"));
	if (flags.contains("use_db_genome")){
		seqgen.useLocalFiles(false);
	}
	else{
		String genomePath = Args.parseString(args, "genome", null);
		if (genomePath!=null){
			seqgen.setGenomePath(genomePath);
		}
	}
	
    // load coordinates
	ArrayList<Region> regions = new ArrayList<Region>();
	if (window==-1)
		regions = CommonUtils.loadRegionFile(Args.parseString(args, "regions", null), genome);
	else{
		ArrayList<Point> points = CommonUtils.loadCgsPointFile(Args.parseString(args, "coords", null), genome);
		for (Point p:points){
			regions.add(p.expand(window/2));
		}
	}
    
    int count=1;
    StringBuilder sb = new StringBuilder();
    for (Region r:regions){
    	String seq = seqgen.execute(r);
    	sb.append(">seq_").append(count).append(" ").append(r.toString()).append("\n");
    	sb.append(seq).append("\n");
    	count++;
    }
    CommonUtils.writeFile(Args.parseString(args, "out", "noname")+"_"+window+"bp.fasta", sb.toString());
    
  }
  
}
