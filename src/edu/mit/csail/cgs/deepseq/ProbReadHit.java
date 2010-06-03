package edu.mit.csail.cgs.deepseq;

import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * A probabilistic Read. This class could be used in cases where a discrete extension of the read is inappropriate.
 * The probabilistic read is implemented by a per-base series of weights over the region surrounding the read.  
 *
 * Not yet implemented!!!
 * 
 * @author shaun
 *
 */
public class ProbReadHit extends ReadHit{

	private BindingModel model;
	
	public ProbReadHit(Genome g, int id,String c, int s, int e, char str, double w, BindingModel bm){
		super(g,id,c,s,e,str,w);
		model=bm;
	}
}
