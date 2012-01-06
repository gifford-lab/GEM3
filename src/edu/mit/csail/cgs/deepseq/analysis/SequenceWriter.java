package edu.mit.csail.cgs.deepseq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.Set;

import net.sf.samtools.util.SequenceUtil;

import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.motifs.WeightMatrix;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.analysis.MultiTF_Binding.Site;
import edu.mit.csail.cgs.deepseq.utilities.CommonUtils;
import edu.mit.csail.cgs.ewok.verbs.SequenceGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;

public class SequenceWriter {

	Genome genome=null;
	private SequenceGenerator<Region> seqgen;
	String[] args;
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		SequenceWriter sw = new SequenceWriter(args);
		sw.writeSequences();

	}

	public SequenceWriter(String[] args){
		this.args = args;
		try {
	    	Pair<Organism, Genome> pair = Args.parseGenome(args);
	    	if(pair==null){
	    	  System.err.println("No genome provided; provide a Gifford lab DB genome name");
	    	  System.exit(1);
	    	}else{
	    		genome = pair.cdr();
	    	}
	    } catch (NotFoundException e) {
	      e.printStackTrace();
	    }

		Set<String> flags = Args.parseFlags(args);
		seqgen = new SequenceGenerator<Region>();
		seqgen.useCache(!flags.contains("no_cache"));
	}
	
	public void writeSequences(){
		int range = Args.parseInteger(args, "range", 101);
		String coorFile = Args.parseString(args, "coords", null);
		ArrayList<String> coorStrs = CommonUtils.readTextFile(coorFile);
		String indexFile = Args.parseString(args, "index", null);
		ArrayList<String> idxStrs = CommonUtils.readTextFile(indexFile);
		StringBuilder sb = new StringBuilder();
		for (String idx: idxStrs){
			if (idx.length()==0)
				continue;
			int i = Integer.parseInt(idx);
			String[] fs = coorStrs.get(i).split("\t");
			if (fs.length==0)
				continue;
			Region r = Point.fromString(genome, fs[0]).expand(range/2);
			String seq = "AACCTTTG";
//			seq = seqgen.execute(r);
			boolean isMinus = fs[1].equals("-1");
			if (isMinus)
				seq = SequenceUtils.reverseComplement(seq);
			sb.append(">").append(r.toString()).append("\t").append(isMinus?"-":"+").append("\n");
			sb.append(seq).append("\n");			
		}
		CommonUtils.writeFile(coorFile+".fasta.txt", sb.toString());
	}
}
