package edu.mit.csail.cgs.deepseq.utilities;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.features.EnrichedFeature;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class SeqEnrichmentByRegion {

	protected DeepSeqExpt signal;
	protected DeepSeqExpt control;
	protected int regionExtension = 1000;
	protected Genome gen=null;
	protected boolean scanGenesOnly=true;
	protected boolean annotOverlapOnly=true;
	protected int maxAnnotDistance=50000;
	protected ArrayList<AnnotationLoader> geneAnnotations = new ArrayList<AnnotationLoader>();
	protected boolean noControl=true;
	protected boolean dbconnected=true;
	protected ArrayList<SeqEnrichResult> results=new ArrayList<SeqEnrichResult>();
	
	public SeqEnrichmentByRegion(String[] args){
		
		ArgParser ap = new ArgParser(args);
		try {
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			gen = pair.cdr();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		List<ChipSeqLocator> dbexpts = Args.parseChipSeq(args,"dbexpt");
        List<ChipSeqLocator> dbctrls = Args.parseChipSeq(args,"dbctrl");
        List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt");
        List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args,"rdbctrl");
        List<File> expts = Args.parseFileHandles(args, "expt");
        List<File> ctrls = Args.parseFileHandles(args, "ctrl");
        boolean nonUnique = ap.hasKey("nonunique") ? true : false;
        String fileFormat = Args.parseString(args, "format", "ELAND");
    	int readLength = Args.parseInteger(args,"readlen",32);
        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
        	signal = new DeepSeqExpt(gen, expts, nonUnique, fileFormat, readLength);
        }else if(dbexpts.size()>0 && expts.size() == 0){
        	signal = new DeepSeqExpt(gen, dbexpts, "db", readLength);
        	dbconnected=true;
        }else if(rdbexpts.size()>0 && expts.size() == 0){
            	signal = new DeepSeqExpt(gen, rdbexpts, "readdb", readLength);
            	dbconnected=true;
        }else{System.err.println("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");System.exit(1);}
        
        if(ctrls.size()>0 && dbctrls.size() == 0 && dbctrls.size()==0){
        	control = new DeepSeqExpt(gen, ctrls, nonUnique, fileFormat, readLength); noControl=false;
        }else if(dbctrls.size()>0 && ctrls.size() == 0){
        	control = new DeepSeqExpt(gen, dbctrls, "db", readLength); noControl=false;
        	dbconnected=true;
        }else if(rdbctrls.size()>0 && ctrls.size() == 0){
        	control = new DeepSeqExpt(gen, rdbctrls, "readdb", readLength); noControl=false;
        	dbconnected=true;
        }else{
        	if(dbctrls.size()>0 && ctrls.size()>0){
        		System.err.println("Cannot mix files and db loading yet...");;System.exit(1);
        	}else{
        		noControl=true; control=null;
        	}
        }
       
	
		//Print some info
		System.out.println("Signal hit count: "+(int)signal.getHitCount()+", weight: "+(int)signal.getWeightTotal());
		if(!noControl)
			System.out.println("Control hit count: "+(int)control.getHitCount()+", weight: "+(int)control.getWeightTotal());
		
        //Gene Annotations
        Collection<String> tfiles = Args.parseStrings(args,"transcripts");
        Collection<String> dbgenes = Args.parseStrings(args,"dbgenes");
        //Special case default
        if(dbgenes.size()==0 && tfiles.size()==0 && dbconnected){
        	String geneSource = gen.getSpecies().equals("Saccharomyces cerevisiae") ? "sgdGene":"refGene";
        	geneAnnotations.add(new AnnotationLoader(gen, geneSource,geneSource, maxAnnotDistance, annotOverlapOnly));
		}
        for(String s:dbgenes)
        	geneAnnotations.add(new AnnotationLoader(gen, s, "refGene", maxAnnotDistance, annotOverlapOnly));
        for(String s:tfiles)
        	geneAnnotations.add(new AnnotationLoader(gen, s, "file", maxAnnotDistance, annotOverlapOnly));		
	}
	
	public void execute(){
		//Naive scaling factor for now 
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		
		Iterator<Gene> testRegions=null;
		//Initialize regions
		ArrayList<Gene> geneList = new ArrayList<Gene>();
		for(AnnotationLoader loader : geneAnnotations){
			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while(chroms.hasNext()){
				NamedRegion c = chroms.next();
				for(Gene r : loader.getGenes(c)){
					geneList.add(r);
                }
			}
		}
		testRegions = geneList.iterator();
		
		while (testRegions.hasNext()) {
			Gene currG = testRegions.next();
			Region currReg = currG.expand(regionExtension, regionExtension);
			List<ReadHit> hits = signal.loadHits(currReg);
			double sigweight = 0, ctrlweight=0;
			for(ReadHit h : hits){	sigweight+=h.getWeight();}
			if(!noControl){
				List<ReadHit> chits = control.loadHits(currReg);
				for(ReadHit h : chits){	ctrlweight+=h.getWeight();}
				ctrlweight*=control.getScalingFactor();
			}
			double sigavg = sigweight>0 ? sigweight/(double)currReg.getWidth() : 0;
			double ctrlavg = ctrlweight>0 ? ctrlweight/(double)currReg.getWidth() : 0;
			double fold = ctrlweight>0 ? sigweight/ctrlweight : sigweight;
			double logfold = fold>0 ? Math.log(fold) : 0;
			double pval = ctrlweight>0 || sigweight>0 ? binomialPValue(ctrlweight, ctrlweight+sigweight) : 0.5;
			if(fold<1)
				pval=1-pval;
			
			SeqEnrichResult r = new SeqEnrichResult();
			r.g=currG;
			r.sig=sigweight;
			r.ctrl=ctrlweight;
			r.fold=fold;
			r.logfold=logfold;
			r.p=pval;
			results.add(r);
		}
		Collections.sort(results);
		results = benjaminiHochbergCorrection(results);
	}
	
	//Binomial CDF assuming scaled control. Uses COLT binomial test
	// k=scaled control, n=scaled control+signal
	protected double binomialPValue(double k, double n){
		double pval=1;
		Binomial b = new Binomial((int)Math.ceil(n), 0.5, new DRand());
		pval = b.cdf((int) Math.ceil(k));
		return(pval);		
	}
	//Multiple hypothesis testing correction -- assumes peaks ordered according to p-value
	protected ArrayList<SeqEnrichResult> benjaminiHochbergCorrection(ArrayList<SeqEnrichResult> x){
		double total = x.size();
		//double significanceThres=0.01;
		ArrayList<SeqEnrichResult> res = new ArrayList<SeqEnrichResult>();
		double rank =1;
		for(SeqEnrichResult y : x){
			y.p = y.p*(total/rank);
			if(y.p>1)
				y.p=1;
			//if(p.score<=significanceThres)
				res.add(y);
			rank++;
		}
		Collections.sort(res);
		return(res);
	}
	protected void printResults(){
		if(!noControl){
			System.out.println("Name\tLength\tSigWeight\tScaledCtrlWeight\tFoldChange\tLog2FoldChange\tP-Value");
		}else{
			System.out.println("Name\tLength\tSigWeight");
		}
		
		for(SeqEnrichResult r : results){
			if(!noControl)
				System.out.println(r.g.getName()+"\t"+r.g.getWidth()+"\t"+r.sig+"\t"+r.ctrl+"\t"+r.fold+"\t"+r.logfold+"\t"+r.p);
			else
				System.out.println(r.g.getName()+"\t"+r.g.getWidth()+"\t"+r.sig);
		}
	}
	
	public static void main(String[] args){
		SeqEnrichmentByRegion sebr = new SeqEnrichmentByRegion(args);
		sebr.execute();
		sebr.printResults();
	}
	
	public class SeqEnrichResult implements Comparable<SeqEnrichResult>{
		public Gene g;
		public double sig, ctrl, fold, logfold, p;		
		
		public int compareTo(SeqEnrichResult x) {
			if(p<x.p){return(-1);}
			else if(p>x.p){return(1);}
			else{return(0);}
		}
	}
}
