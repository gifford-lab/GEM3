package edu.mit.csail.cgs.deepseq.discovery;

import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.utils.Pair;


public abstract class MultiConditionFeatureFinder extends FeatureFinder{
	
	/**
	 * It is a <tt>List</tt> containing the experiments used. <br>
	 * Each element of the <tt>List</tt> represents a condition of an IP-Control 
	 * experiment <tt>Pair</tt>. <br>
	 * So, the size of the <tt>List</tt> gives the number of conditions. <br>
	 * Each <tt>Pair</tt> of IP-Control experiment contains the ChipSeq data for the
	 * specified condition. Each experiment is stored in the form of <tt>DeepSeqExpt</tt>
	 * and can hold multiple replicates of it (through the <tt>ReadLoader</tt> field).
	 */
	protected ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> experiments = new ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>>();
	
	protected int numConditions=0;
	
	
	public MultiConditionFeatureFinder(Genome g, ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> expts) {
		super(g);
		numConditions=expts.size();
	}
	public MultiConditionFeatureFinder(String[] args){
		super(args);
	}
	public MultiConditionFeatureFinder(String[] args, Genome g, ArrayList<Pair<DeepSeqExpt,DeepSeqExpt>> expts){
		super(args);
		gen = g;
		genomeLen = g.getGenomeLength();
		experiments = expts;
	}

	//Call this method before exiting
	public void cleanup(){
		for(Pair<DeepSeqExpt,DeepSeqExpt> p : experiments){
			if(p.car()!=null)
				p.car().closeLoaders();
			if(p.cdr()!=null)
				p.cdr().closeLoaders();
		}
	}
}
