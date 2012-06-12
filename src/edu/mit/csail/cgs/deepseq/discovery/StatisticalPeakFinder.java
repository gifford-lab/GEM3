package edu.mit.csail.cgs.deepseq.discovery;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;
import java.util.Collection;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import cern.jet.stat.Probability;

import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.deepseq.BackgroundCollection;
import edu.mit.csail.cgs.deepseq.BindingModel;
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
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.data.DataFrame;
import edu.mit.csail.cgs.utils.models.data.DataRegression;

/**
 * Superclass for all peak finders that use the binomial test for population proportion equality.
 * 
 * @author shaun
 *
 */
public abstract class StatisticalPeakFinder extends SingleConditionFeatureFinder{

	protected BackgroundCollection signalBacks=new BackgroundCollection(), ctrlBacks=new BackgroundCollection();
	protected BackgroundCollection signalPerBaseBack=new BackgroundCollection(), ctrlPerBaseBack=new BackgroundCollection();
	
    protected double minFoldChange = 1;
	protected double binWidth=100;
	protected double binStep = 25;
	protected double highLogConf=-9; 
	protected double perBaseLogConf=-7;
	protected int fixedPerBaseCutoff = -1; ////Ignore more than this number of reads mapping to a base (-1 turns feature off)
	protected double significanceThres=0.01;
	protected double read5PrimeExt=50, read3PrimeExt=50,readShift=50;
    protected boolean towerfiltering=true, needlefiltering=true;
	protected int towerWindow=(int)binWidth*10;
	protected List<Integer> dbacks=new ArrayList<Integer>(); 
	protected int MAXSECTION = 100000000;
	protected boolean peakLRBal=false, peakWithModel=false; //Place peaks at max count, bindingModels or LR balance point??? 
	protected BindingModel bindingModel=null; //only used for placing peaks
	protected boolean buildBindingModel=false;
	protected List<EnrichedFeature> signalPeaks = new Vector<EnrichedFeature>();
	protected List<EnrichedFeature> controlPeaks = new Vector<EnrichedFeature>();
	protected int scalingWindow = 10000; 
	protected boolean addAllToScaling=false;
	protected boolean useBinomialTest=true;
	protected boolean multiHypoTest=true;
	protected boolean printProgress=true;
	protected boolean showGeneAnnotations=true;
	protected boolean showOtherAnnotations=true;
	protected boolean countReads = true; //use this if you expect the per-base thresholds will make a significant difference

    private Binomial binomial = new Binomial(10,.5, new DRand());
    private int longestRead = 100;
	
	//Constructors
	public StatisticalPeakFinder(DeepSeqExpt signal){this(signal,null);}	
	public StatisticalPeakFinder(DeepSeqExpt signal, DeepSeqExpt ctrl){
		super(signal, ctrl);
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		initializeBackgrounds();	
	}
	public StatisticalPeakFinder(String[] args){
		super(args);
		
		//Initialize scaling factor
		if(!noControl)
			control.setScalingFactor(signal.getWeightTotal()/control.getWeightTotal());
		
		//Read extensions
		read5PrimeExt = Args.parseDouble(args,"read5ext",read5PrimeExt);
		read3PrimeExt = Args.parseDouble(args,"read3ext",read3PrimeExt);
		readShift = Args.parseDouble(args,"readshift",readShift);
				
		//Load options for peak calling 
        setMinFoldChange(Args.parseDouble(args,"min_fold_change",minFoldChange));
        setBinWidth(Args.parseDouble(args,"binwidth",binWidth));
        setBinStep(Args.parseDouble(args,"binstep",binStep));
        setHighLogConf(Args.parseDouble(args,"highlogconf",highLogConf));
        setSigThres(Args.parseDouble(args,"sigthres",significanceThres));
        setPerBaseLogConf(Args.parseDouble(args,"pblogconf",perBaseLogConf));
        setFixedPerBase(Args.parseInteger(args,"fixedpb",fixedPerBaseCutoff));
        setCountReads(Args.parseFlags(args).contains("countreads"));
        setTowerFilter(!Args.parseFlags(args).contains("allowtowers"));
        setNeedleFilter(!Args.parseFlags(args).contains("allowneedles"));
        setShowGeneAnnotations(!Args.parseFlags(args).contains("noShowGeneAnnotations"));
        setShowOtherAnnotations(!Args.parseFlags(args).contains("noShowOtherAnnotations"));        
        MAXSECTION = Args.parseInteger(args,"maxloadlen",MAXSECTION);
       // metaPeak = Args.parseArgs(args).contains("printpeakdistrib");
        dbacks= (List<Integer>) Args.parseIntegers(args, "dynback");
        showGeneAnnotations = !Args.parseFlags(args).contains("nogeneannots");
        showOtherAnnotations = !Args.parseFlags(args).contains("nootherannots");
        if(!stranded){
        	setLRBalPeaks(Args.parseFlags(args).contains("balpeak"));
        	setModelPeaks(Args.parseArgs(args).contains("modelpeak"));
        	loadBindingModel(Args.parseString(args,"model",null));
        	buildModel(Args.parseFlags(args).contains("buildmodel"));
        	if(bindingModel != null){
			    //readShift = (double)bindingModel.maxShift();
			    //read5PrimeExt = (double)bindingModel.maxShift()-(readLength/2);
			    //read3PrimeExt = (double)bindingModel.maxShift()-(readLength/2);
        		readShift = 0;
        		Pair<Double, Double> probInterval = bindingModel.probIntervalAboveUniform();
			    read5PrimeExt = -1*probInterval.car();
			    read3PrimeExt = probInterval.cdr()-readLength;
			}
		}
        
        signal.setFivePrimeExt((int)read5PrimeExt); 
		signal.setThreePrimeExt((int)read3PrimeExt); 
		signal.setShift((int)readShift); 
		if(!noControl){
			control.setFivePrimeExt((int)read5PrimeExt);
			control.setThreePrimeExt((int)read3PrimeExt);
			control.setShift((int)readShift);
		}
	}
	
	//Run the peak caller
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
		signalFeatures.addAll(signalPeaks);
		controlFeatures.addAll(controlPeaks);
		
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
	 * To be counted as a potential peak, the read counts in a window must be over all "signalBackgrounds" and under all "controlBackgrounds".
	 * By default, the background collections contain only a global Poisson model. 
	 * If the use of local backgrounds is requested (including gene only scanning), the following occurs:
	 * 		If there is a control, local low thresholds are added to the control
	 * 		If there is a control and the control is scaled, a local high threshold based on expectation from CONTROL is added to the signal
	 * 		If there is no control or the control is not scaled, a local high threshold based on expectation from SIGNAL is added to the signal
	 * This makes the behavior of local background thresholds similar to that used by MACS. 
	 */
	protected void initializeBackgrounds(){
		signalBacks=new BackgroundCollection(); ctrlBacks=new BackgroundCollection();
		signalPerBaseBack=new BackgroundCollection(); ctrlPerBaseBack=new BackgroundCollection();
		
		//Set up the background models
        double readTotalLen = readLength +read5PrimeExt+read3PrimeExt;
        if(scanGenesOnly){dbacks.add(0);}

        //Non-stranded
        if(!stranded){
            //Per Base
            signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, signal.getWeightTotal(), 1, genomeLen, mappableGenome, 1, 1, '.', 1, true));
            if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, control.getWeightTotal(), 1, genomeLen, mappableGenome, 1, 1, '.', 1, true));}
            
        	//Signal genomic
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, signal.getWeightTotal(), readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '.', 1, true));
        	//Control genomic
    		if(!noControl){
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, control.getWeightTotal(), readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '.', 1, true));
    		}
            //Locals
            for(Integer i : dbacks){
            	if(!noControl){ //local control low 
                	//signal threshold based on what would be expected from the CONTROL locality
                	signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, control.getWeightTotal(), readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '.', control.getScalingFactor(), false));
            		ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, signal.getWeightTotal(), readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '.', 1/control.getScalingFactor(), false));
            	}else{
            		//local signal high -- this may bias against locally enriched signal regions, and so should only be used if there is no control or if the control is not yet scaled
                	if(i.intValue()>=5000) // we don't want the window too small in this case
                		signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, signal.getWeightTotal(), readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '.', 1, true));
            	}	
            }            
        }else{ //Stranded
        	double sigPosW = signal.getStrandedWeightTotal('+');
        	double sigNegW = signal.getStrandedWeightTotal('-');
        	double ctrlPosW = control.getStrandedWeightTotal('+');
        	double ctrlNegW = control.getStrandedWeightTotal('-');
            //Per Base
            signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, sigPosW, 1, genomeLen, mappableGenome, 1, 1, '+', 1, true));
            signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, sigNegW, 1, genomeLen, mappableGenome, 1, 1, '-', 1, true));
            if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, ctrlPosW, 1, genomeLen, mappableGenome, 1, 1, '+', 1, true));}
            if(!noControl){ctrlPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, ctrlNegW, 1, genomeLen, mappableGenome, 1, 1, '-', 1, true));}
            
        	//Signal genomic high & low (low for when signal & control are swapped)
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, sigPosW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '+', 1, true));
        	signalBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, sigNegW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '-', 1, true));
        	//Control high & low (high for when signal & control are swapped)
    		if(!noControl){
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, ctrlPosW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '+', 1, true));
    			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(-1, highLogConf, ctrlNegW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '-', 1, true));
    		}
            //Locals
            for(Integer i : dbacks){
            	if(!noControl){ //local control low 
                	//signal threshold based on what would be expected from the CONTROL locality
            		signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, ctrlPosW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '+', control.getScalingFactor(), false));
            		signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, ctrlNegW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '-', control.getScalingFactor(), false));
            		if(!noControl){
            			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, sigPosW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '+', 1/control.getScalingFactor(), false));
            			ctrlBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, sigNegW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '-', 1/control.getScalingFactor(), false));
            		}
            	}else{
            		//local signal high -- this may bias against locally enriched signal regions, and so should only be used if there is no control or if the control is not yet scaled
            		signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, sigPosW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '+', 1, true));
            		signalBacks.addBackgroundModel(new PoissonBackgroundModel(i.intValue(), highLogConf, sigNegW, readTotalLen, genomeLen, mappableGenome, binWidth, binStep, '-', 1, true));
            	}
            }
        }
	}
	
    class StatisticalThread implements Runnable {
        private Collection<Region> regions;
        //These two structures are used by makeHitLandscape
        private double[] landscape=null;
        private double[] startcounts=null;
        private int numStrandIter;
        private StatisticalPeakFinder parent;
        private double ipTotHits, backTotHits;
        private boolean recordForScaling, postProcess;
        private List<PairedCountData> scalingPairs;

        public StatisticalThread(Collection<Region> r, 
                                 List<PairedCountData> scalingPairs,
                                 int numStrandIter, 
                                 double ipTotHits, 
                                 double backTotHits, 
                                 boolean recordForScaling,
                                 boolean postProcess,
                                 StatisticalPeakFinder parent) {
            regions = r;
            this.scalingPairs = scalingPairs;
            this.numStrandIter = numStrandIter;
            this.ipTotHits = ipTotHits;
            this.backTotHits = backTotHits;
            this.recordForScaling = recordForScaling;
            this.postProcess = postProcess;
            this.parent = parent;
        }
        //Makes integer arrays corresponding to the read landscape over the current region
        protected void makeHitLandscape(ArrayList<ReadHit> hits, Region currReg, int perBaseMax, char strand){
            int numBins = (int)(currReg.getWidth()/binStep);
            int [] counts = new int[currReg.getWidth()+1];
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
        }

        public void run() {
            for (Region currentRegion : regions) {
                double printStep=50000000,  numPrint=1;	   
                double sigHitsScalingWindow=0, ctrlHitsScalingWindow=0;
                boolean peakInScalingWindow=false;
                EnrichedFeature lastSigPeak=null, lastCtrlPeak=null;
                Region currScalingRegion=null;
                int basesDone = 0;
                //Split the job up into large chunks
                for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION){
                    int y = x+MAXSECTION; 
                    if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                    Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
                    ArrayList<ReadHit> ipHits = new ArrayList<ReadHit>();
                    ArrayList<ReadHit> backHits = new ArrayList<ReadHit>();
				
                    synchronized(signal) {
                        ipHits.addAll(signal.loadExtHits(currSubRegion));
                    }
                    if (!noControl) {
                        synchronized(control) {
                            backHits.addAll(control.loadExtHits(currSubRegion));
                        }
                    }

                    for(int stranditer=1; stranditer<=numStrandIter; stranditer++){
                        System.err.println("Working on " + currSubRegion + ", " + stranditer);
                        ArrayList<EnrichedFeature> currSigRes = new ArrayList<EnrichedFeature>();
                        ArrayList<EnrichedFeature> currCtrlRes = new ArrayList<EnrichedFeature>();
                        lastSigPeak=null; lastCtrlPeak=null;
                        //If stranded peak-finding, run over both strands separately
                        char str = !stranded ? '.' : (stranditer==1 ? '+' : '-');
					
                        int signalPB = fixedPerBaseCutoff>0 ? fixedPerBaseCutoff : signalPerBaseBack.getMaxThreshold(str); 
                        makeHitLandscape(ipHits, currSubRegion, signalPB, str);
                        double ipStackedHitCounts[] = landscape.clone();
                        double ipHitStartCounts[] = startcounts.clone();
                        double backStackedHitCounts[] = null;
                        double backHitStartCounts[] = null;
                        if (!noControl) {
                            makeHitLandscape(backHits, currSubRegion, ctrlPerBaseBack.getMaxThreshold(str), str);
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
                            if(signalBacks.passesGenomicThreshold((int)ipWinHits, str) && (noControl || ctrlBacks.underGenomicThreshold((int)backWinHits, str))){
                                //Second Test: refresh all thresholds & test again
                                signalBacks.updateModels(currSubRegion, i-x, ipStackedHitCounts, backStackedHitCounts);
                                if(!noControl)
                                    ctrlBacks.updateModels(currSubRegion, i-x, backStackedHitCounts, ipStackedHitCounts);
                                //System.out.println("UPDATED:"+i);
                                //signalBacks.printThresholds();
                                if(signalBacks.passesAllThresholds((int)ipWinHits, str) && (noControl || ctrlBacks.underAllThresholds((int)backWinHits, str))){
                                    Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+binWidth-1));
                                    //Add hit to signalPeaks
                                    lastSigPeak=addEnrichedReg(currSigRes, lastSigPeak, currWin, ipWinHits, (noControl ? 0 : control.getScalingFactor()*backWinHits), ipTotHits, (noControl ? 0 : control.getScalingFactor()*backTotHits), str);
                                    peakInScalingWindow=true;
                                } else{
                                    sigHitsScalingWindow+=ipHitStartCounts[currBin];
                                }
                            } else {
                                sigHitsScalingWindow+=ipHitStartCounts[currBin];
                            }

						
                            //Ctrl vs Signal Enrichment 
                            if(!noControl){
                                //First Test: is the read count above the genome-wide thresholds? 
                                if(ctrlBacks.passesGenomicThreshold((int)backWinHits, str) && (signalBacks.underGenomicThreshold((int)ipWinHits, str))){
                                    //Second Test: refresh all thresholds & test again
                                    signalBacks.updateModels(currSubRegion, i-x, ipStackedHitCounts, backStackedHitCounts);
                                    ctrlBacks.updateModels(currSubRegion, i-x, backStackedHitCounts, ipStackedHitCounts);
                                    if(ctrlBacks.passesAllThresholds((int)backWinHits, str) && (signalBacks.underAllThresholds((int)ipWinHits, str))){
                                        Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+binWidth-1));
                                        //Add hit to controlPeaks
                                        lastCtrlPeak=addEnrichedReg(currCtrlRes, lastCtrlPeak, currWin, control.getScalingFactor()*backWinHits, ipWinHits, control.getScalingFactor()*backTotHits, ipTotHits, str);
                                        peakInScalingWindow=true;
                                    }else{
                                        ctrlHitsScalingWindow+=backHitStartCounts[currBin];
                                    }
                                }else{
                                    ctrlHitsScalingWindow+=backHitStartCounts[currBin];
                                }
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
                            if(printProgress){
                                basesDone+=binStep;
                                if(basesDone >= numPrint*printStep){
                                    System.out.println(String.format("%s, %d bases, strand %d", currentRegion.toString(), basesDone, stranditer));
                                                                     
                                    numPrint++;
                                }
                            }
                            currBin++;
                        }
                        //Now count the total reads for each peak
                        for (ReadHit h : ipHits) {
                            int w = h.getWidth();
                            if (w > longestRead) {
                                longestRead = w;
                            }
                        }
                        Collections.sort(ipHits);
                        for (ReadHit h : backHits) {
                            int w = h.getWidth();
                            if (w > longestRead) {
                                longestRead = w;
                            }
                        }
                        Collections.sort(backHits);

                        countTotalReadsInPeaks(currSigRes, ipHits, backHits, true);
                        countTotalReadsInPeaks(currCtrlRes, backHits, ipHits, false);
					
                        if(postProcess){
                            //Trim                            
                            trimPeaks(currSigRes, ipHits,str);
                            trimPeaks(currCtrlRes, backHits,str);

                            //filter towers?
                            if(towerfiltering){
                                ArrayList<EnrichedFeature> tmpSigRes = filterTowers(currSigRes, currCtrlRes);
                                ArrayList<EnrichedFeature> tmpCtrlRes = filterTowers(currCtrlRes, currSigRes);
                                currSigRes = tmpSigRes; currCtrlRes = tmpCtrlRes; 
                            }
						
                            //Find the peaks
                            for(EnrichedFeature h : currSigRes){
                                if(peakLRBal)
                                    h.peak = findPeakLRBalance(ipHits, h.coords);
                                else if(peakWithModel && bindingModel!=null)
                                    h.peak = findPeakWithBindingModel(ipHits, h.coords, bindingModel);
                                else {
                                    h.peak = findPeakMaxHit(ipHits, h.coords,str);
                                }
                                

                            }
                        }
                        //Add to results
                        signalPeaks.addAll(currSigRes);
                        controlPeaks.addAll(currCtrlRes);
					
                    }// end of for(int stranditer=1; stranditer<=numStrandIter; stranditer++) loop
				
                }// end of for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION) loop
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
		//Peak placing info
		if(peakLRBal)
			System.out.println("Peaks will be placed by L-R balancing");
		else if(peakWithModel){
			if(bindingModel==null)
				System.out.println("Peaks should be placed by binding model, but no binding model provided");
			else
				System.out.println("Peaks will be placed by binding model");
		}else
			System.out.println("Peaks will be placed by maximum overlap");
		
		signalPeaks = new Vector<EnrichedFeature>();
		controlPeaks = new Vector<EnrichedFeature>();
		List<PairedCountData> scalingPairs = new Vector<PairedCountData>();

		double ipTotHits = countReads ? reCountReads(signal, signalPerBaseBack) : signal.getWeightTotal();
        double backTotHits = noControl ? 1 : control.getWeightTotal(); // set to 1 if no background so the binomial test doesn't return NaN        
		int numStrandIter = stranded ? 2 : 1; 

        Thread[] threads = new Thread[maxThreads];
        ArrayList threadRegions[] = new ArrayList[maxThreads];
        int i = 0;
        for (i = 0 ; i < threads.length; i++) {
            threadRegions[i] = new ArrayList<Region>();
        }
        while (testRegions.hasNext()) {
            threadRegions[(i++) % maxThreads].add(testRegions.next());
        }

        for (i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new StatisticalThread(threadRegions[i], 
                                                        scalingPairs,
                                                        numStrandIter, 
                                                        ipTotHits,
                                                        backTotHits,
                                                        recordForScaling,
                                                        postProcess,
                                                        this));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) { }
            for (i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }
        }
        		
		//Sort & correct for multiple hypotheses
		Collections.sort(signalPeaks);
		Collections.sort(controlPeaks);
		if(multiHypoTest){
			signalPeaks=benjaminiHochbergCorrection(signalPeaks);
			Collections.sort(signalPeaks);
			controlPeaks=benjaminiHochbergCorrection(controlPeaks);
			Collections.sort(controlPeaks);
		}
		
		if(recordForScaling && !noControl){
			double s = estimateScalingFactor(scalingPairs);
			System.out.println(String.format("Scaling factor = %.3f",s));
			control.setScalingFactor(s);
		}		
	}

	
	//Accessors
    public void setMinFoldChange(double f) {minFoldChange = f;}
	public void setBinWidth(double w){binWidth = w;}
	public void setBinStep(double s){binStep = s;}
	public void setReadLength(int l){readLength=l;}
	public void setReadShift(int s){readShift=(double)s;signal.setShift(s); if(!noControl){control.setShift(s);}}
	public void set5PrimeExt(int e){read5PrimeExt=(double)e;signal.setFivePrimeExt(e); if(!noControl){control.setFivePrimeExt(e);}}
	public void set3PrimeExt(int e){read3PrimeExt=(double)e;signal.setThreePrimeExt(e); if(!noControl){control.setThreePrimeExt(e);}}
	public void setHighLogConf(double p){highLogConf=p;}
	public void setPerBaseLogConf(double p){perBaseLogConf=p;}
	public void setFixedPerBase(int f){fixedPerBaseCutoff=f;}
	public void setSigThres(double s){significanceThres=s;}
	public void setTowerFilter(boolean f){towerfiltering=f;}
	public void setCountReads(boolean cr){countReads = cr;}
	public void setNeedleFilter(boolean f){needlefiltering=f;}
	public void setLRBalPeaks(boolean lr){peakLRBal=lr;}
	public void setModelPeaks(boolean m){peakWithModel=m;}
	public void buildModel(boolean b){buildBindingModel=b;}
	public void setShowGeneAnnotations(boolean ga){showGeneAnnotations=ga;}
	public void setShowOtherAnnotations(boolean oa){showOtherAnnotations=oa;}
	public void loadBindingModel(String m){if(m!=null){bindingModel = new BindingModel(new File(m));}}
	public void setBindingModel(BindingModel bm){
		if(bm != null){
			bindingModel = bm;
			readShift = 0;
    		Pair<Double, Double> probInterval = bindingModel.probIntervalAboveUniform();
		    read5PrimeExt = -1*probInterval.car();
		    read3PrimeExt = probInterval.cdr()-readLength;
		    signal.setFivePrimeExt((int)read5PrimeExt); 
			signal.setThreePrimeExt((int)read3PrimeExt); 
			signal.setShift((int)readShift); 
			if(!noControl){
				control.setFivePrimeExt((int)read5PrimeExt);
				control.setThreePrimeExt((int)read3PrimeExt);
				control.setShift((int)readShift);
			}
		}
	}
	//"Manually" count reads 
	protected double reCountReads(DeepSeqExpt e, BackgroundCollection back){
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
	
	//count the total reads within each peak.  ipHits and backHits must be sorted first
	protected void countTotalReadsInPeaks(ArrayList<EnrichedFeature> peaks, ArrayList<ReadHit> ipHits, ArrayList<ReadHit> backHits, boolean forIPPeaks){
		for(EnrichedFeature peak : peaks){
			peak.signalTotalHits= forIPPeaks ? overlappingHits(ipHits, peak.coords, peak.strand).size() : overlappingHits(ipHits, peak.coords, peak.strand).size()*control.getScalingFactor();
			peak.backTotalHits = (backHits==null || backHits.size()==0) ? 0 : (forIPPeaks ? overlappingHits(backHits, peak.coords, peak.strand).size()*control.getScalingFactor() : overlappingHits(backHits, peak.coords, peak.strand).size());
			peak.overrep = peak.backTotalHits==0 ? peak.signalTotalHits : peak.signalTotalHits/peak.backTotalHits;
		}
	}
	
    /* hits must be sorted and longestRead must be set from it */
	protected ArrayList<ReadHit> overlappingHits(ArrayList<ReadHit> hits, Region window, char str){
		ArrayList<ReadHit> sub = new ArrayList<ReadHit>();
        int l = 0;
        int r = hits.size();
        int windowStart = window.getStart();
        int windowEnd = window.getEnd();
        while (r - l > 10) {
            int c = (l + r) / 2;
            if (windowStart > hits.get(c).getStart()) {
                l = c;
            } else {
                r = c;
            }
        }
        while (l > 0 && (hits.get(l).getStart() + longestRead > windowStart)) {
            l--;
        }
        while (l < hits.size() && hits.get(l).getStart() <= windowEnd) {
            ReadHit hit = hits.get(l);
			if(window.overlaps(hit) && (str=='.' || str==hit.getStrand())){sub.add(hit);}			
            l++;
        }
		return(sub);
	}
	
	//Add a peak to the pile
	protected EnrichedFeature addEnrichedReg(ArrayList<EnrichedFeature> currres, EnrichedFeature lastHit, Region currWin, double ipWinHits, double backWinHits, double ipTotHits, double backTotHits, char str){
		EnrichedFeature resHit=null;
		double score = -1;
		if(useBinomialTest)
			score = binomialPValue(backWinHits, ipWinHits+backWinHits);
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
	// k=scaled control, n=scaled control+signal
	protected double binomialPValue(double k, double n){
		double pval=1;
        synchronized (binomial) {
            binomial.setNandP((int)Math.ceil(n), 1.0 / (minFoldChange + 1));
            pval = binomial.cdf((int) Math.ceil(k));
        }
		return(pval);		
	}
	//Multiple hypothesis testing correction -- assumes peaks ordered according to p-value
	protected ArrayList<EnrichedFeature> benjaminiHochbergCorrection(List<EnrichedFeature> peaks){
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
		}
        return(res);
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
	
	
	// Filter the towers; protected because this method only makes sense from within the callEnrichedRegions method 
	protected ArrayList<EnrichedFeature> filterTowers(ArrayList<EnrichedFeature> thisChannel, ArrayList<EnrichedFeature> otherChannel){
		ArrayList<EnrichedFeature> res = new ArrayList<EnrichedFeature>();
		double length = 2*(readLength); 
		for(EnrichedFeature peak : thisChannel){
			boolean isTower=false; //innocent until proven guilty
			
			//Filter 1: remove "peaks" that neighbor background towers
			//Given a local background, I'd be surprised if this comes into play at all
			for(EnrichedFeature tower : otherChannel){
				if(peak.coords.getChrom().equals(tower.coords.getChrom()) && peak.coords.distance(tower.coords)<=towerWindow)
					isTower=true;
			}
			//Filter 2: simple minimum size filter. Peaks should be more than two (unextended) reads long, right?!?
			if(peak.coords.getWidth()<=length)
				isTower=true;
			
			
			if(!isTower){
				res.add(peak);
			}
		}
		return(res);
	}	
	
	//Estimate the scaling factor 
	protected double estimateScalingFactor(List<PairedCountData> scalingData){
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
	/* Trim the peaks back to the coordinates of the first & last read in the peak. 
       ipHits must be sorted
     */
	protected void trimPeaks(ArrayList<EnrichedFeature> curr, ArrayList<ReadHit> ipHits, char str){
		for(EnrichedFeature peak : curr){
			ArrayList<ReadHit> winHits = overlappingHits(ipHits, peak.coords,str);
			if(winHits.size()>0){
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
	}
	/* Find the exact peak locations based on maximum overlapping read counts. 
       Reads must be sorted and longestRead set *before* calling this method
       Careful; make sure the reads are extended... */
	protected Point findPeakMaxHit(ArrayList<ReadHit> hits, Region coords, char str){
		ArrayList<ReadHit> winHits = overlappingHits(hits, coords,str);
		int [] sum = new int [coords.getWidth()+1];
		for(int s=0; s<sum.length; s++){sum[s]=0;}

		for(ReadHit r : winHits){
            int start = r.getStart()-coords.getStart(); 
            int stop= r.getEnd()-coords.getStart();
            for(int i=start; i<stop; i++){
                if(i>=0 && i<sum.length){
                    sum[i]++;
                }
            }
		}
		int max = 0, maxPos = -1;
		for(int s=0; s<sum.length; s++){
			if(sum[s]>max){
				max= sum[s];
				maxPos=s;
			}
		}
		Point p = new Point(gen, coords.getChrom(), maxPos+coords.getStart());
		return(p);
	}
	/* Find the exact peak locations based on left-right balance of forward/reverse reads */
	protected Point findPeakLRBalance(ArrayList<ReadHit> hits, Region coords){
		int [] forward = new int [coords.getWidth()+1];
		int [] reverse = new int [coords.getWidth()+1];
		for(int s=0; s<=coords.getWidth(); s++){forward[s]=0; reverse[s]=0;}
		int ext = (int)read5PrimeExt;
		for(StrandedRegion r : hits){
			int readStart = r.getStrand()=='+' ? r.getFivePrime()+ext : r.getFivePrime()-ext;
			if(readStart>=coords.getStart() && readStart<=coords.getEnd()){
				if(r.getStrand()=='+'){
					forward[Math.max(0, r.getStart()-coords.getStart())]++;
				}else{
					reverse[Math.min(r.getEnd()-coords.getStart(),coords.getWidth())]++;
				}
			}
		}
		int minBal = 10000, minPos = -1;
		for(int s=0; s<=coords.getWidth(); s++){
			int left=0, right=0;
			for(int l=0; l<=s; l++){left+=forward[l];}
			for(int r=coords.getWidth(); r>s; r--){right+=reverse[r];}
			if(Math.abs(left-right)<minBal){
				minBal =Math.abs(left-right);
				minPos = s;
			}
		}
		Point p = new Point(gen, coords.getChrom(), minPos+coords.getStart());
		return(p);
	}
	/* Find exact peak using a BindingModel */
	protected Point findPeakWithBindingModel(ArrayList<ReadHit> hits, Region coords, BindingModel model){
		int maxPos=0; double maxScore=0;
		double [] p = new double[coords.getWidth()+1];
		for(int k=0; k<=coords.getWidth(); k++){p[k]=0;}
		int ext = (int)read5PrimeExt;
		for(ReadHit x : hits){
			int readStart = x.getStrand()=='+' ? x.getFivePrime()+ext : x.getFivePrime()-ext;
			if(readStart>=coords.getStart()-model.getMax() && readStart<=coords.getEnd()+model.getMax()){
				int offset = readStart-coords.getStart();
				if(x.getStrand()=='+')
					for(int i=Math.max(model.getMin()+offset, 0); i<=Math.min(coords.getWidth(), offset+model.getMax()); i++)
						p[i]+=model.probability(i-offset);
				else
					for(int i=Math.max(offset-model.getMax(), 0); i<=Math.min(coords.getWidth(), offset-model.getMin()); i++)
						p[i]+=model.probability(offset-i);
			}
		}
		for(int k=0; k<=coords.getWidth(); k++)
			if(p[k]>maxScore){maxScore=p[k]; maxPos=k;}
			
		Point pt = new Point(gen, coords.getChrom(), coords.getStart()+maxPos);
		return(pt);
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
                "  --format <ELAND/NOVO/BOWTIE/BED (default ELAND)> \n" +
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
                "  --minfold <minimumfoldchange> \n"+           
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
