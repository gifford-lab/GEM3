package edu.mit.csail.cgs.deepseq.discovery;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import cern.jet.stat.Probability;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.deepseq.BackgroundCollection;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.PairedCountData;
import edu.mit.csail.cgs.deepseq.PoissonBackgroundModel;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.features.ClipSeqPeak;
import edu.mit.csail.cgs.deepseq.features.EnrichedFeature;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.ewok.verbs.ChromRegionIterator;
import edu.mit.csail.cgs.ewok.verbs.ChromosomeGenerator;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.models.data.DataRegression;

/**
 * Superclass for all peak finders that use the binomial test for population proportion equality.
 * 
 * @author shaun
 *
 */
public abstract class DifferentialEnrichmentFinder extends SingleConditionFeatureFinder{

	protected BackgroundCollection signalBacks=new BackgroundCollection(), ctrlBacks=new BackgroundCollection();
	//protected BackgroundCollection signalPerBaseBack=new BackgroundCollection(), ctrlPerBaseBack=new BackgroundCollection();
	
	protected double binWidth=100;
	protected double binStep = 25;
	protected double highLogConf=-9; 
	protected double perBaseLogConf=-7;
	protected int fixedPerBaseCutoff = -1; ////Ignore more than this number of reads mapping to a base (-1 turns feature off)
	protected double significanceThres=0.01;
	protected double foldDiffThreshold = 1.5;
	protected int MAXSECTION = 100000000;
	protected ArrayList<EnrichedFeature> signalEvents = new ArrayList<EnrichedFeature>();
	protected ArrayList<EnrichedFeature> controlEvents = new ArrayList<EnrichedFeature>();
	protected int scalingWindow = 10000; 
	protected boolean needlefiltering=false;
	protected boolean addAllToScaling=false;
	protected boolean useBinomialTest=true;
	protected boolean multiHypoTest=true;
	protected boolean printProgress=true;
	protected boolean showGeneAnnotations=true;
	protected boolean showOtherAnnotations=true;
	private double[] landscape=null;
	private double[] startcounts=null;
	
	//Constructors
	public DifferentialEnrichmentFinder(DeepSeqExpt signal, DeepSeqExpt ctrl){
		super(signal, ctrl);
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		initializeBackgrounds();	
	}
	public DifferentialEnrichmentFinder(String[] args){
		super(args);
		
		//Initialize scaling factor
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		
		//Load options for differential event calling 
        setBinWidth(Args.parseDouble(args,"binwidth",binWidth));
        setBinStep(Args.parseDouble(args,"binstep",binStep));
        setHighLogConf(Args.parseDouble(args,"highlogconf",highLogConf));
        setSigThres(Args.parseDouble(args,"sigthres",significanceThres));
        setPerBaseLogConf(Args.parseDouble(args,"pblogconf",perBaseLogConf));
        setFixedPerBase(Args.parseInteger(args,"fixedpb",fixedPerBaseCutoff));
        setShowGeneAnnotations(!Args.parseFlags(args).contains("noShowGeneAnnotations"));
        setShowOtherAnnotations(!Args.parseFlags(args).contains("noShowOtherAnnotations"));        
        MAXSECTION = Args.parseInteger(args,"maxloadlen",MAXSECTION);

	}
	
	//Run the differential event caller
	//This method may be over-ridden by specific peak-finding implementations
	public List<Feature> execute(){
		
		//Initialize background models
		initializeBackgrounds();
		
		//Initialize the iterator for the test regions/genes
		Iterator<Region> testRegions=null;
		if(!scanGenesOnly)
			testRegions = new ChromosomeGenerator().execute(gen);
		else{
			ArrayList<Region> geneList = new ArrayList<Region>();
			for(AnnotationLoader loader : geneAnnotations){
				ChromRegionIterator chroms = new ChromRegionIterator(gen);
				while(chroms.hasNext()){
					NamedRegion c = chroms.next();
					for(Region r : loader.getAnnotations(c)){
						if(r.getWidth()>=(binWidth*2))
							geneList.add(r);
                    }
				}
			}
			testRegions = geneList.iterator();
		}
		
		//Call the peak finder
		callEnrichedRegions(testRegions, true, true);
		signalFeatures.addAll(signalEvents);
		controlFeatures.addAll(controlEvents);
		
		//Add closest genes
		System.out.println("Adding gene annotations");
		addClosestGenes(signalFeatures);
		
		//Add other annotations
		System.out.println("Adding other annotations");
		addRegionAnnotations(signalFeatures);
		
		return signalFeatures;
	}
	
	public void orderFeatsByLocation() {
		if(signalFeatures.size() > 0)
			signalFeatures = orderByLocation(signalFeatures);
		
		if(controlFeatures.size() > 0)
			controlFeatures = orderByLocation(controlFeatures);
	}//end of orderFeatsByLocation method
	
	private List<Feature> orderByLocation(List<Feature> list) {
		
		Map<String, Feature> map = new TreeMap<String, Feature>();
		List<Feature> orderFeats = new ArrayList<Feature>();
		
		for(Feature feat:list)
			map.put(feat.coords.toString(), feat);
		
		for(String key:map.keySet())
			orderFeats.add(map.get(key));
		
		return orderFeats;
	}//end of orderByLocation method
	
	
	/**
	 * Set up the background models
	 * 
	 * Here's the background policy:
	 * To be counted as a potential events in the signal, the read counts in a window must be over all "signalBackgrounds".
	 * By default, the background collections contain only a global Poisson model. 
	 */
	protected void initializeBackgrounds(){
		signalBacks=new BackgroundCollection(); ctrlBacks=new BackgroundCollection();
		//signalPerBaseBack=new BackgroundCollection(); ctrlPerBaseBack=new BackgroundCollection();
		
		
        //Non-stranded
        if(!stranded){
            //Per Base
           //signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, signal.getWeightTotal(), 1, genomeLen, mappableGenome, 1, 1, '.', 1, true));
           // if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, control.getWeightTotal(), 1, genomeLen, mappableGenome, 1, 1, '.', 1, true));}
            
        	//Signal genomic
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, signal.getWeightTotal(), readLength, genomeLen, mappableGenome, binWidth, binStep, '.', 1, true));
        	//Control genomic
    		if(!noControl){
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, control.getWeightTotal(), readLength, genomeLen, mappableGenome, binWidth, binStep, '.', 1, true));
    		}
                   
        }else{ //Stranded
        	double sigPosW = signal.getStrandedWeightTotal('+');
        	double sigNegW = signal.getStrandedWeightTotal('-');
        	double ctrlPosW = control.getStrandedWeightTotal('+');
        	double ctrlNegW = control.getStrandedWeightTotal('-');
            //Per Base
           // signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, sigPosW, 1, genomeLen, mappableGenome, 1, 1, '+', 1, true));
           // signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, sigNegW, 1, genomeLen, mappableGenome, 1, 1, '-', 1, true));
           // if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, ctrlPosW, 1, genomeLen, mappableGenome, 1, 1, '+', 1, true));}
           // if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, ctrlNegW, 1, genomeLen, mappableGenome, 1, 1, '-', 1, true));}
            
        	//Signal genomic high & low (low for when signal & control are swapped)
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, sigPosW, readLength, genomeLen, mappableGenome, binWidth, binStep, '+', 1, true));
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, sigNegW, readLength, genomeLen, mappableGenome, binWidth, binStep, '-', 1, true));
        	//Control high & low (high for when signal & control are swapped)
    		if(!noControl){
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, ctrlPosW, readLength, genomeLen, mappableGenome, binWidth, binStep, '+', 1, true));
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, ctrlNegW, readLength, genomeLen, mappableGenome, binWidth, binStep, '-', 1, true));
    		}
        }
	}
	
	/**
	 * The main functional implementation for the StatisticalPeakFinder
	 * Returns a list of Peaks passing the threshold
	 * @param testRegions
	 * @return
	 */
	public void callEnrichedRegions(Iterator<Region> testRegions, boolean postProcess, boolean recordForScaling){
		double basesDone=0, printStep=10000000,  numPrint=0;;
		
		signalEvents = new ArrayList<EnrichedFeature>();
		controlEvents = new ArrayList<EnrichedFeature>();
		EnrichedFeature lastSigPeak=null, lastCtrlPeak=null;
		ArrayList<PairedCountData> scalingPairs=new ArrayList<PairedCountData>();
		Region currScalingRegion=null;
		boolean peakInScalingWindow=false;
		double sigHitsScalingWindow=0, ctrlHitsScalingWindow=0;
		double ipTotHits = signal.getWeightTotal();
        double backTotHits = noControl ? 1 : control.getWeightTotal(); // set to 1 if no background so the binomial test doesn't return NaN        
		int numStrandIter = stranded ? 2 : 1; 
		if(printProgress){System.out.print("Progress (bp): ");}
		
		while (testRegions.hasNext()) {
			Region currentRegion = testRegions.next();
			
			lastSigPeak=null; lastCtrlPeak=null;
			//Split the job up into large chunks
			for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION){
				int y = x+MAXSECTION; 
				if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
				Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
				ArrayList<ReadHit> ipHits = new ArrayList<ReadHit>();
				ArrayList<ReadHit> backHits = new ArrayList<ReadHit>();
				
				ipHits.addAll(signal.loadExtHits(currSubRegion));
				if (!noControl)
					backHits.addAll(control.loadExtHits(currSubRegion));

				for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
					ArrayList<EnrichedFeature> currSigRes = new ArrayList<EnrichedFeature>();
					ArrayList<EnrichedFeature> currCtrlRes = new ArrayList<EnrichedFeature>();
					lastSigPeak=null; lastCtrlPeak=null;
					//If stranded peak-finding, run over both strands separately
					char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					
					int signalPB = 1;//fixedPerBaseCutoff>0 ? fixedPerBaseCutoff : signalPerBaseBack.getMaxThreshold(str); 
					makeHitLandscape(ipHits, currSubRegion, signalPB, str);
					double ipStackedHitCounts[] = landscape.clone();
					double ipHitStartCounts[] = startcounts.clone();
					double backStackedHitCounts[] = null;
					double backHitStartCounts[] = null;
	                if (!noControl) {
	                	//makeHitLandscape(backHits, currSubRegion, ctrlPerBaseBack.getMaxThreshold(str), str);
	                	makeHitLandscape(backHits, currSubRegion, 1, str);
	                	backStackedHitCounts = landscape.clone();
	            		backHitStartCounts = startcounts.clone();
	                }
					
	                //initialize scaling stuff
	                currScalingRegion = new Region(gen, currSubRegion.getChrom(), currSubRegion.getStart(), currSubRegion.getStart()+scalingWindow<=currSubRegion.getEnd() ? currSubRegion.getStart()+scalingWindow-1 : currSubRegion.getEnd());
	                sigHitsScalingWindow=0; ctrlHitsScalingWindow=0;
	                peakInScalingWindow=false;
	                //Scan regions
	                int currBin=0;
					for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)binWidth; i+=(int)binStep){
						double ipWinHits=ipStackedHitCounts[currBin];
						double backWinHits= noControl ? 0 : backStackedHitCounts[currBin];
						
						//Signal vs Ctrl Enrichment 
						//First Test: is the read count above the genome-wide thresholds? 
						if(signalBacks.passesGenomicThreshold((int)ipWinHits, str) && ipWinHits/(control.getScalingFactor()*backWinHits) >=  foldDiffThreshold){
							Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+binWidth-1));
							//Add hit to signalPeaks
							lastSigPeak=addEnrichedReg(currSigRes, lastSigPeak, currWin, ipWinHits, (noControl ? 0 : backWinHits), ipTotHits, (noControl ? 0 : backTotHits), str);
							peakInScalingWindow=true;
						}else{sigHitsScalingWindow+=ipHitStartCounts[currBin];}
						
						//Ctrl vs Signal Enrichment 
						if(!noControl){
						//First Test: is the read count above the genome-wide thresholds? 
						if(ctrlBacks.passesGenomicThreshold((int)backWinHits, str) && (control.getScalingFactor()*backWinHits)/ipWinHits >=  foldDiffThreshold){
							Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+binWidth-1));
							//Add hit to controlPeaks
							lastCtrlPeak=addEnrichedReg(currCtrlRes, lastCtrlPeak, currWin, backWinHits, ipWinHits, backTotHits, ipTotHits, str);
							peakInScalingWindow=true;
						}else{ctrlHitsScalingWindow+=backHitStartCounts[currBin];}
						}
								
						//WARNING: In this loop, the peak's total read counts and over-representation will not be accurate 

						//Scaling
						if(recordForScaling && !noControl){
							if(i>currScalingRegion.getEnd()){//Gone past end of scaling window... reset
								if(currScalingRegion.getWidth()==scalingWindow){
									if(addAllToScaling || !peakInScalingWindow)
										scalingPairs.add(new PairedCountData(sigHitsScalingWindow, ctrlHitsScalingWindow));
								}
								peakInScalingWindow=false;
								sigHitsScalingWindow=0;
								ctrlHitsScalingWindow=0;
								currScalingRegion = new Region(gen, currSubRegion.getChrom(), i, i+scalingWindow<=currSubRegion.getEnd() ? i+scalingWindow-1 : currSubRegion.getEnd());
							}
						}
						
						//Print out progress
						if(stranditer==1 && printProgress){
							basesDone+=binStep;
							if(basesDone > numPrint*printStep){
								if(numPrint%10==0){System.out.print(String.format("(%.0f)", (numPrint*printStep)));}
								else{System.out.print(".");}
								if(numPrint%50==0 && numPrint!=0){System.out.print("\n");}
								numPrint++;
							}
						}
						currBin++;
					}
					//Now count the total reads for each peak
					countTotalReadsInPeaks(currSigRes, ipHits, backHits, true);
					countTotalReadsInPeaks(currCtrlRes, backHits, ipHits, false);
					
					if(postProcess){
						//Trim
						trimPeaks(currSigRes, ipHits,str);
						trimPeaks(currCtrlRes, backHits,str);						
					}
					//Add to results
					signalEvents.addAll(currSigRes);
					controlEvents.addAll(currCtrlRes);
					
				}// end of for(int stranditer=1; stranditer<=numStrandIter; stranditer++) loop
				
			}// end of for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION) loop
			
		}// end of while (testRegions.hasNext()) loop
		
		if(printProgress){System.out.println(String.format("(%.0f)", basesDone));}
		//Sort & correct for multiple hypotheses
		Collections.sort(signalEvents);
		Collections.sort(controlEvents);
		if(multiHypoTest){
			signalEvents=benjaminiHochbergCorrection(signalEvents);
			Collections.sort(signalEvents);
			controlEvents=benjaminiHochbergCorrection(controlEvents);
			Collections.sort(controlEvents);
		}
		
		if(recordForScaling && !noControl){
			double s = estimateScalingFactor(scalingPairs);
			System.out.println(String.format("Scaling factor = %.3f",s));
			control.setScalingFactor(s);
		}
		
	}

	
	//Accessors
	public void setBinWidth(double w){binWidth = w;}
	public void setBinStep(double s){binStep = s;}
	public void setReadLength(int l){readLength=l;}
	public void setHighLogConf(double p){highLogConf=p;}
	public void setPerBaseLogConf(double p){perBaseLogConf=p;}
	public void setFixedPerBase(int f){fixedPerBaseCutoff=f;}
	public void setSigThres(double s){significanceThres=s;}
	public void setShowGeneAnnotations(boolean ga){showGeneAnnotations=ga;}
	public void setShowOtherAnnotations(boolean oa){showOtherAnnotations=oa;}
	
	//"Manually" count reads 
/*	protected double reCountReads(DeepSeqExpt e, BackgroundCollection back){
		double count =0;
		int perBaseMax=0;
		ChromRegionIterator chroms = new ChromRegionIterator(gen);
		while(chroms.hasNext()){
			NamedRegion c = chroms.next();
			for(int x=c.getStart(); x<=c.getEnd(); x+=MAXSECTION){
				int y = x+MAXSECTION; 
				if(y>c.getEnd()){y=c.getEnd();}
				Region currSubRegion = new Region(gen, c.getChrom(), x, y);
				
				char str = '+';
				perBaseMax = fixedPerBaseCutoff>0 ? fixedPerBaseCutoff : signalPerBaseBack.getMaxThreshold(str);
				for(Float f : signal.loadStrandedBaseCounts(currSubRegion, str).cdr()){
					if(f<=perBaseMax){count+=f;}
					else{count+=perBaseMax;}
				}
				str = '-';
				perBaseMax = fixedPerBaseCutoff>0 ? fixedPerBaseCutoff : signalPerBaseBack.getMaxThreshold(str);
				signal.loadStrandedBaseCounts(currSubRegion, str);
				for(Float f : signal.loadStrandedBaseCounts(currSubRegion, str).cdr()){
					if(f<=perBaseMax){count+=f;}
					else{count+=perBaseMax;}
				}
			}
		}
		System.err.println("Recounted experiment: total weight = "+count);
		return(count);
	}
*/
	
	//Makes integer arrays corresponding to the read landscape over the current region
	protected void makeHitLandscape(ArrayList<ReadHit> hits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth()/binStep);
		int [] counts = new int[currReg.getWidth()+1];
		//double [] land = new double[numBins+1];
		landscape = new double[numBins+1];
		startcounts = new double[numBins+1];
		for(int i=0; i<=numBins; i++){landscape[i]=0; startcounts[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(ReadHit r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getStart()-currReg.getStart(),0,currReg.getWidth());
				counts[offset]++;//small issue here... counts will not be corrected for scalingFactor in control
				if(!needlefiltering || (counts[offset] <= perBaseMax)){
					int binstart = inBounds((int)((double)offset/binStep), 0, numBins);
					int binend = inBounds((int)((double)(r.getEnd()-currReg.getStart())/binStep), 0, numBins);
					for(int i=binstart; i<=binend; i++){
						landscape[i]+=r.getWeight();
					}
					if(r.getStrand()=='+')
						startcounts[binstart]+=r.getWeight();
					else
						startcounts[binend]+=r.getWeight();
				}
			}
		}
		//return(land);
	}
	
	//count the total reads within each peak 
	protected void countTotalReadsInPeaks(ArrayList<EnrichedFeature> peaks, ArrayList<ReadHit> ipHits, ArrayList<ReadHit> backHits, boolean forIPPeaks){
		for(EnrichedFeature peak : peaks){
			peak.signalTotalHits= forIPPeaks ? overlappingHits(ipHits, peak.coords, peak.strand).size() : overlappingHits(ipHits, peak.coords, peak.strand).size()*control.getScalingFactor();
			peak.backTotalHits = (backHits==null || backHits.size()==0) ? 0 : (forIPPeaks ? overlappingHits(backHits, peak.coords, peak.strand).size()*control.getScalingFactor() : overlappingHits(backHits, peak.coords, peak.strand).size());
			peak.overrep = peak.backTotalHits==0 ? peak.signalTotalHits : peak.signalTotalHits/peak.backTotalHits;
		}
	}
	
	protected ArrayList<ReadHit> overlappingHits(ArrayList<ReadHit> hits, Region window, char str){
		ArrayList<ReadHit> sub = new ArrayList<ReadHit>();
		for(ReadHit r : hits){
			if(window.overlaps(r) && (str=='.' || str==r.getStrand())){sub.add(r);}			
		}
		return(sub);
	}
	
	//Add a peak to the pile
	protected EnrichedFeature addEnrichedReg(ArrayList<EnrichedFeature> currres, EnrichedFeature lastHit, Region currWin, double ipWinHits, double backWinHits, double ipTotHits, double backTotHits, char str){
		EnrichedFeature resHit=null;
		double score = -1;
		if(useBinomialTest)
			score = binomialPValue(ipWinHits, ipTotHits, backWinHits/backTotHits);
		else
			score = binomialSampleEquality(ipWinHits, backWinHits, ipTotHits, backTotHits);
		
		//Is this hit close to the previously added one?
		if(lastHit!=null && currWin.distance(lastHit.coords)<=binWidth){
			lastHit.coords=new Region(gen, lastHit.coords.getChrom(), lastHit.coords.getStart(), currWin.getEnd()-1);
			
			if(score<lastHit.score || (noControl && ipWinHits>lastHit.signalMaxHits)){
				lastHit.score=score;
				lastHit.signalMaxHits=ipWinHits;
				lastHit.backMaxHits=backWinHits;
				lastHit.overrep = backWinHits > 0 ? (ipWinHits/ipTotHits)/(backWinHits/backTotHits) :-1;
			}
			resHit=lastHit;
		}else{
			double over = backWinHits > 0 ? (ipWinHits/ipTotHits)/(backWinHits/backTotHits) :-1;
			if(stranded){
				ClipSeqPeak hit = new ClipSeqPeak(currWin, ipWinHits, backWinHits, score, over, str);
				currres.add(hit);
				resHit=hit;
			}else{
				EnrichedFeature hit = new EnrichedFeature(currWin, ipWinHits, backWinHits, score, over);
				currres.add(hit);
				resHit=hit;											
			}										
		}return(resHit);
	}
	//Binomial CDF assuming scaled control. Uses COLT binomial test
	// k=signal in region, n=total signal, p = scaled control/scaled total control
	protected double binomialPValue(double k, double n, double p){
		double pval=1;
		Binomial b = new Binomial((int)Math.ceil(n), p, new DRand());
		pval = 1 - b.cdf((int) Math.ceil(k));
		return(pval);		
	}
	//Multiple hypothesis testing correction -- assumes peaks ordered according to p-value
	protected ArrayList<EnrichedFeature> benjaminiHochbergCorrection(ArrayList<EnrichedFeature> peaks){
		double total = peaks.size();
		ArrayList<EnrichedFeature> res = new ArrayList<EnrichedFeature>();
		double rank =1;
		for(EnrichedFeature p : peaks){
			p.score = p.score*(total/rank);
			if(p.score>1)
				p.score=1;
			if(p.score<=significanceThres)
				res.add(p);
			rank++;
		}return(res);
	}
	// Binomial test for differences between two population proportions 
	protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		if(!Double.isNaN(Z))
			return(1-Probability.normal(Z));
		else
			return(-1);
	}
	
	
	//Estimate the scaling factor 
	protected double estimateScalingFactor(ArrayList<PairedCountData> scalingData){
		double slope = 0;
		if(scalingData==null || scalingData.size()==0)
			return 1;
		DataFrame df = new DataFrame(edu.mit.csail.cgs.deepseq.PairedCountData.class, scalingData.iterator());
		DataRegression r = new DataRegression(df, "x~y - 1");
		r.calculate();
		Map<String, Double> map = r.collectCoefficients();
		slope = map.get("y");
		return(slope);
	}
	/* Trim the peaks back to the coordinates of the first & last read in the peak*/
	protected void trimPeaks(ArrayList<EnrichedFeature> curr, ArrayList<ReadHit> ipHits, char str){
		for(EnrichedFeature peak : curr){
			ArrayList<ReadHit> winHits = overlappingHits(ipHits, peak.coords,str);
			StrandedRegion min=winHits.get(0);
			StrandedRegion max=winHits.get(0);
			for(StrandedRegion sr : winHits){
				if(sr.getStart()<min.getStart()){min=sr;}
				if(sr.getEnd()>max.getEnd()){max=sr;}
			}
			int startOff = peak.coords.getStart()-min.getStart()<0 ? peak.coords.getStart()-min.getStart(): 0;
			int endOff = max.getEnd()-peak.coords.getEnd()<0 ? max.getEnd()-peak.coords.getEnd():0;

			peak.coords = peak.coords.expand(startOff, endOff);
		}
	}

	//keep the number in bounds
	protected final double inBounds(double x, double min, double max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	public void printArgs(){
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
                "  --format <ELAND/NOVO/BOWTIE/BED/SAM/TOPSAM (default ELAND)> \n" +
                "  --nonunique [use nonunique reads]\n" +
                "Required:\n"+
                "  --species <organism name;genome version>\n  OR\n"+
                "  --geninfo <file with chr name/length pairs> \n" +
                "Options:\n"+
                "  --mappable <mappable proportion of the genome (default 0.8) \n" +
                "  --out <output file root name> \n"+
                "  --seqwin <sequence length> \n" +
                "  --noseqs [flag to not print sequences] \n"+
                "  --binwidth <width of bins> \n"+
                "  --binstep <offset of bins> \n"+
                "  --sigthres <significance threshold \n"+
                "  --highlogconf <log 10 Poisson threshold for signal channel> \n" +
                "  --lowlogconf <log 10 Poisson threshold for control channel> \n" +
                "  --pblogconf <log 10 Poisson threshold per base> \n" +
                "  --fixedpb <fixed per base occurrence> \n" +
                "  --dynback <dynamic background thresholds (-1/0/window = genomic/gene/local window)> \n" +
                "  --readlen <length> \n" +
                "  --read5ext <5' extension> --read3ext <3' extension> --readshift <shift tags this distance> \n"+
                "  --balpeak [place peaks at forward-reverse balance point] \n" +
                "  --modelpeak [place peaks with model (model reqd)] \n" +
                "  --model <binding model file>\n" +
                "  --buildmodel [re-build the binding model] \n" +
                "  --allowtowers [don't filter towers] \n" +
                "  --allowneedles [don't filter needles] \n" +
                "  --stranded [search each strand separately] \n" +
                "  --dbgenes <name of gene annotation to examine in Gifford DB> \n" +
                "  --printgff [print GFF format output\n" +
                "Other DB annotation sources: --namedregions --namedstrandedregions --namedtypedregions --repeatmasker\n" +
                "  --transcripts <flat-file with gene annotations> \n" +
                "  --annot <flat-file with other annotations>\n" +
                "  --maxannotdist <distance from gene> \n" +
                "  --annotoverlap [only look at annotations that overlap a peak] \n" +
                "  --scangenesonly [only look for peaks within genes (CLIP-seq specific)] \n" +
                "  --maxloadlen <load reads from this length region each iteration; set low if running out of memory>\n");
	}
	
}
