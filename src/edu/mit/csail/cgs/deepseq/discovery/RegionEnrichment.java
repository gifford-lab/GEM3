package edu.mit.csail.cgs.deepseq.discovery;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.features.EnrichedNamedRegion;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.tools.utils.Args;

public class RegionEnrichment extends SingleConditionFeatureFinder {
	public ArrayList<EnrichedNamedRegion> geneMeasurements;
	public int regionExtension = 1000;
	public boolean scanPromotersOnly=false;
	
	//Constructors
	public RegionEnrichment(DeepSeqExpt signal){this(signal,null);}	
	public RegionEnrichment(DeepSeqExpt signal, DeepSeqExpt ctrl){
		super(signal, ctrl);
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
	}
	public RegionEnrichment(String[] args){
		super(args);
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		scanOnlyPromoters(Args.parseFlags(args).contains("scanpromotersonly"));
		setRegionExt(Args.parseInteger(args, "regionext", 1000));
	}
	
	public List<Feature> execute() {
		//Initialize the iterator for the test regions/genes
		Iterator<NamedRegion> testRegions=null;
		if(!scanGenesOnly && !scanPromotersOnly){
			//testRegions = new ChromosomeGenerator().execute(gen);
			//CHANGE TO LOAD REGIONS FROM FILE (remember to name the regions)
		}else{
			ArrayList<NamedRegion> geneList = new ArrayList<NamedRegion>();
			for(AnnotationLoader loader : geneAnnotations){
				ChromRegionIterator chroms = new ChromRegionIterator(gen);
				while(chroms.hasNext()){
					NamedRegion c = chroms.next();
					for(Gene r : loader.getGenes(c)){
						int start=-1, end=-1;
						if(scanGenesOnly){
							start = r.getStart()-regionExtension < 1 ? 1 : r.getStart()-regionExtension;
							end = r.getEnd()+regionExtension > c.getEnd() ? c.getEnd() : r.getEnd()+regionExtension;
						}else if(scanPromotersOnly){
							start = r.getTSS()-regionExtension < 1 ? 1 : r.getTSS()-regionExtension;
							end = r.getTSS()+regionExtension > c.getEnd() ? c.getEnd() : r.getTSS()+regionExtension;
						}
						geneList.add(new NamedRegion(r.getGenome(), r.getChrom(), start, end, r.getName()));
                    }
				}
			}
			testRegions = geneList.iterator();
		}
		
		//Call the enrichment analyzer
		analyzeEnrichment(testRegions);

		signalFeatures.addAll(geneMeasurements);

		return(signalFeatures);
	}
	
	
	public void analyzeEnrichment(Iterator<NamedRegion> testRegions){
		geneMeasurements = new ArrayList<EnrichedNamedRegion>();
		
		double ipTot = signal.getWeightTotal();
		double ctrlTot = noControl ? 1 : control.getWeightTotal(); // set to 1 if no background so the binomial test doesn't return NaN        
		
		while (testRegions.hasNext()) {
			NamedRegion currReg = testRegions.next();
		
			double ipW = signal.sumWeights(currReg);
			double ctrlW = noControl ? 0.0 : control.sumWeights(currReg)*control.getScalingFactor();
			
			double over = ctrlW>0 ? ipW/ctrlW : -1;
			double ipFPKM = ((ipW/ipTot)*1000000)/((double)currReg.getWidth()/1000);
			double ctrlFPKM = noControl ? 0 : (ctrlW/(ctrlTot*control.getScalingFactor())*1000000)/((double)currReg.getWidth()/1000);
			
			double pval = (ipW+ctrlW)>0 ? binomialPValue(ctrlW, ipW+ctrlW) : 1.0;
			
			EnrichedNamedRegion currRes = new EnrichedNamedRegion(currReg, ipW, ctrlW, ipFPKM, ctrlFPKM, pval, over);
			
			//Add to results
			geneMeasurements.add(currRes);
			
		}// end of while (testRegions.hasNext()) loop
		
		
		//Sort & correct for multiple hypotheses
		Collections.sort(geneMeasurements);
		geneMeasurements=benjaminiHochbergCorrection(geneMeasurements);
		Collections.sort(geneMeasurements);
		
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
	protected ArrayList<EnrichedNamedRegion> benjaminiHochbergCorrection(ArrayList<EnrichedNamedRegion> x){
		double total = x.size();
		ArrayList<EnrichedNamedRegion> res = new ArrayList<EnrichedNamedRegion>();
		double rank =1;
		for(EnrichedNamedRegion y : x){
			y.score = y.score*(total/rank);
			if(y.score>1)
				y.score=1;
			res.add(y);
			rank++;
		}
		Collections.sort(res);
		return(res);
	}
	
	public static void main(String[] args){
		RegionEnrichment re = new RegionEnrichment(args);
		List<Feature> res = re.execute();
		re.printFeatures();
	}
	
	//Accessors
	public void scanOnlyPromoters(boolean s){scanPromotersOnly=s;}
	public void setRegionExt(int e){regionExtension=e;}

	public void printError() {
		System.err.println("RegionEnrichment: analyze ChIP-seq read counts over entire genes or regions.\n");
		System.err.println("Usage:\n " +
				"Using with Gifford Lab ReadDB:\n" +
                "  --rdbexpt <solexa expt> \n" +
                "  --rdbctrl <background expt> \n" +
                "Using with Gifford Lab Oracle DB:\n" +
                "  --dbexpt <solexa expt> \n" +
                "  --dbctrl <background expt> \n" +
                "Using with flat-files:\n" +
                "  --expt <aligned reads file for expt> \n" +
                "  --ctrl <aligned reads file for ctrl> \n" +
                "  --format <ELAND/NOVO/BOWTIE/BED (default ELAND)> \n" +
                "  --nonunique [use nonunique reads]\n" +
                "Required:\n"+
                "  --species <organism name;genome version>\n  OR\n"+
                "  --geninfo <file with chr name/length pairs> \n" +
                "Options:\n" +
                "  --scanpromotersonly [scan RefGene promoters]\n" +
                "  --scangenesonly [scan RefGene gene regions]\n" +
                "  --scanregions <file of coordinates> NOT YET IMPLEMENTED!\n" +
                "  --regionext <bp to extend BOTH sides>" +
                "");
	}
	
}
