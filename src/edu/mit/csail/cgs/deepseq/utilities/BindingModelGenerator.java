package edu.mit.csail.cgs.deepseq.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.mit.csail.cgs.datasets.chipseq.ChipSeqLocator;
import edu.mit.csail.cgs.datasets.general.Point;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.BackgroundModel;
import edu.mit.csail.cgs.deepseq.BindingModel;
import edu.mit.csail.cgs.deepseq.DeepSeqExpt;
import edu.mit.csail.cgs.deepseq.PoissonBackgroundModel;
import edu.mit.csail.cgs.deepseq.ReadHit;
import edu.mit.csail.cgs.deepseq.features.EnrichedFeature;
import edu.mit.csail.cgs.deepseq.features.EnrichedFeatureFileReader;
import edu.mit.csail.cgs.ewok.verbs.PointParser;
import edu.mit.csail.cgs.ewok.verbs.RegionParser;
import edu.mit.csail.cgs.projects.readdb.PairedHit;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

public class BindingModelGenerator {
	private BindingModel model;
	private Genome gen =null;
	private List<EnrichedFeature> peaks=null;
	private int window=1000;
	private int min=-2000, max=2000;
	private double overrepFilter=5;
	private double pFilter=0.01;
	private DeepSeqExpt IP=null;
	private List<EnrichedFeature> towers = new ArrayList<EnrichedFeature>();
	private boolean smooth=true; 
	private int smoothWin=5;//10bp smoothing
	private int binSize=5; //size of bin to take 
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") ||(!ap.hasKey("expt")&&!ap.hasKey("rdbexpt")) || !ap.hasKey("peaks")) { 
            System.err.println("Usage:\n " +
                   "BindingModelGenerator \n" +
                   " Required: \n" +
                   "  --species <species;genome> \n" +
                   "  --rdbexpt <IP expt names> \n" +
                   "  OR \n" +
                   "  --expt <alignment file> AND --format <ELAND/NOVO/BOWTIE>\n" +
                   "  --paired [flag for paired-end data]\n" +
                   "  --nonunique [use the non-unique hits]\n" +
                   "  --peaks <file containing coordinates of peaks in Shaun's format> \n" +
                   "  --out <output file name>\n" +
                   "  --win <window around peaks to take>\n" +
                   "  --bin <bin size>\n" +
                   "  --splinesmooth [smooth the probabilities with a cubic spline]\n" +
                   "  --p <p-value threshold on input peaks>\n" +
                   "  --or <over-rep threshold on input peaks>\n" +
                   "  --towers <file containing tower locations in Shaun's format>\n" +
                   "");
            return;
        }
        String peaksFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks"):null;
    	int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():1000;
    	int bin = ap.hasKey("bin") ? new Integer(ap.getKeyValue("bin")).intValue():5;
    	double p = ap.hasKey("p") ? new Double(ap.getKeyValue("p")):0.0001;
    	double or = ap.hasKey("or") ? new Double(ap.getKeyValue("or")):5;
        List<ChipSeqLocator> dbexpts = Args.parseChipSeq(args,"dbexpt");
        List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt");
        List<File> expts = Args.parseFileHandles(args, "expt");
        String fileFormat = Args.parseString(args, "format", "ELAND");
        String outName = Args.parseString(args, "out", "out.model");
        boolean nonUnique = ap.hasKey("nonunique") ? true : false;
        boolean splineSmooth = ap.hasKey("splinesmooth") ? true : false;
        String towerFile = ap.hasKey("towers") ? ap.getKeyValue("towers") : null;
    	int rL = Args.parseInteger(args,"readlen",32);
    	boolean pairedEnd = ap.hasKey("paired"); 
    	
        try {
        	Pair<Organism, Genome> pair = Args.parseGenome(args);
			DeepSeqExpt signal=null;
	        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
	        	signal = new DeepSeqExpt(pair.cdr(), expts, nonUnique, fileFormat, rL);
	        }else if(dbexpts.size()>0 && expts.size() == 0){
	        	signal = new DeepSeqExpt(pair.cdr(), dbexpts, "db", rL);	        	
	        }else if(rdbexpts.size()>0 && expts.size() == 0){
	        	signal = new DeepSeqExpt(pair.cdr(), rdbexpts, "readdb", rL);	        	
	        }else{System.err.println("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");System.exit(1);}
        	signal.setPairedEnd(pairedEnd);
        	
			//initialize
			BindingModelGenerator generator = new BindingModelGenerator(pair.cdr(), signal);
			generator.setBinSize(bin);
			generator.setWindow(win);
			generator.setPFilter(p);
			generator.setOverRepFilter(or);
			generator.initPeaks(peaksFile);
			if(towerFile != null)
				generator.loadTowers(towerFile);
			
			//execute
			generator.execute();
			
			if(splineSmooth)
				generator.smoothModel();
			
			//Print the meta peak
			generator.printModel(outName);
			
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public BindingModelGenerator(Genome g, DeepSeqExpt sig){
		IP=sig;
		gen=g;
	}
	public BindingModelGenerator(Genome g, DeepSeqExpt sig, List<EnrichedFeature> p){this(g, sig, p, null);}
	public BindingModelGenerator(Genome g, DeepSeqExpt sig, List<EnrichedFeature> p, List<EnrichedFeature> t){
		IP=sig;
		gen=g;
		peaks = p;
		towers=t;
		System.out.println("BindingModel generating from "+peaks.size()+" peaks.");
	}

	public void initPeaks(String fname){
		peaks = loadPeaks(fname);
		System.out.println("BindingModel generating from "+peaks.size()+" peaks.");
	}
	//Accessor
	public void setWindow(int w){window=w;}
	public void setBinSize(int b){binSize=b;}
	public void setPFilter(double p){pFilter=p;}
	public void setOverRepFilter(double o){overrepFilter=o;}
	public BindingModel getBindingModel(){return model;}
	public void smoothModel(){model.smooth(BindingModel.SMOOTHING_STEPSIZE, BindingModel.SMOOTHING_AVG_PTS);}
	
	//Execute the binding model generator
	//Includes tower & needle filter
	public BindingModel execute(){
		List<Pair<Integer, Double>> bindingDist = new ArrayList<Pair<Integer,Double>>();
		int numBins = ((max-min)/binSize);
		BackgroundModel back = new PoissonBackgroundModel(-1, Math.pow(10, -9), IP.getWeightTotal(), 1, gen.getGenomeLength(), 0.8, 1, 1, '.', 1, true);
		int perBaseMax = back.getThreshold();
		int winLen=max-min;
		double [] forward = new double [numBins+1];
		double [] reverse = new double [numBins+1];
		for(int w=0; w<=numBins; w++){
			forward[w]=0; reverse[w]=0;
		}
		
		System.out.println("Binding model gen filters: p<="+pFilter+" overRep>="+overrepFilter);
		//Load the hits
		int numProc=0;
		int peaksUsed=0, readsUsed=0;
		for(EnrichedFeature p : peaks){
			if(p.score<=pFilter && (p.overrep==-1 || p.overrep >=overrepFilter)){
				peaksUsed++;
				List<ReadHit> hits =null;
				if(!IP.isPairedEnd())
					hits = IP.loadHits(p.coords);
				else{
					hits = new ArrayList<ReadHit>();
					List<PairedHit> phits = IP.loadPairsAsPairs(p.coords);
					HashMap<PairedHit, Double> pairFilter = new HashMap<PairedHit,Double>();
					//Filter for valid pairs
					for(PairedHit ph : phits){
						//Hard filter for any pairs that share both coordinates
						if(!pairFilter.containsKey(ph)){
							if(gen.getChromName(ph.leftChrom).equals(p.coords.getChrom()) && gen.getChromName(ph.rightChrom).equals(p.coords.getChrom()) && 
									!ph.lesserStrand() && ph.greaterStrand()){
								hits.add(convertToReadHit(gen, -1, ph, true));
								hits.add(convertToReadHit(gen, -1, ph, false));
							}
						}
						pairFilter.put(ph, new Double(ph.weight));
					}
				}
				ArrayList<Region> currTowers = new ArrayList<Region>();
				for(EnrichedFeature t : towers)
					if(p.coords.overlaps(t.coords))
						currTowers.add(t.coords);
				int [] counts = new int[p.coords.getWidth()+1];
				for(int i=0; i<=p.coords.getWidth(); i++){counts[i]=0;}
				
				for(ReadHit h : hits){
					readsUsed++;
					boolean inTower=false;
					for(Region t : currTowers)
						if(p.coords.overlaps(t))
							inTower=true;
					if(!inTower){
						int offset=h.getFivePrime()-p.coords.getStart()<0 ? 0 : h.getFivePrime()-p.coords.getStart();
						if(offset>p.coords.getWidth())
							offset=p.coords.getWidth();
						counts[offset]++;
						if(counts[offset] <= perBaseMax){
							int dist = h.getStrand()=='+' ? h.getFivePrime()-p.peak.getLocation() : p.peak.getLocation()-h.getFivePrime();
							if(dist>=min && dist<=max){
								if(h.getStrand()=='+')
									forward[(dist-min)/binSize]++;
								else
									reverse[(dist-min)/binSize]++;
							}
						}
					}
				}
			}
			numProc++;
			//if(numProc%100==0)
			//	System.out.print("\n");
		}//System.out.print("\n");
		System.out.println(peaksUsed+" peaks and "+readsUsed+" reads used to generate binding model");
		
		//Print the distribs & find maxes
		double maxFor=0, maxRev=0;
		double maxForPos=0, maxRevPos=winLen;
		//System.out.println("Tag distributions around peaks");
		//System.out.println("RelPos\tForward\tReverse");
		for(int w=0; w<=numBins; w++){
			int rel = min+(w*binSize);
			
			//System.out.println(rel+"\t"+forward[w]+"\t"+reverse[w]);
			if(forward[w]>maxFor && rel<-10){maxFor=forward[w]; maxForPos=rel;}
			if(reverse[w]>maxRev && rel>10){maxRev=reverse[w];maxRevPos=rel;}		
		}double diff = (maxRevPos-maxForPos)/2;
		System.out.println("Estimated Forward-Reverse Shift:\t"+(int)diff);		
		
		for(int w=0; w<=numBins; w++){
			int rel = min+(w*binSize);
			if(smooth){
				double tot=0; double num=0;
				for(int v=w-(smoothWin/binSize); v<=w+(smoothWin/binSize); v++){
					if(v>=0 && v<winLen){
						tot+=(forward[v]+reverse[v]);
						num++;
					}
				}bindingDist.add(new Pair<Integer,Double>(rel, tot/num));
			}else{
				bindingDist.add(new Pair<Integer,Double>(rel, (forward[w]+reverse[w])));
			}
		}
		model = new BindingModel(bindingDist);
		return model;
	}	
	
	//Load a set of regions from a peak file
	public ArrayList<EnrichedFeature> loadPeaks(String filename){
		ArrayList<EnrichedFeature> pset = new ArrayList<EnrichedFeature>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = reader.readLine();//ignore first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            	            
	            if(words.length>=8 && window!=-1){	//Shaun's format
	            	RegionParser parser = new RegionParser(gen);
	            	Region reg = parser.execute(words[0]);
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
	            	double over = new Double(words[9]);
	            	double pval = new Double(words[6]);
	            	
	            	if(pval<=pFilter && (over==-1 || over >=overrepFilter)){
		            	int rstart = p.getLocation()-(window/2)<1 ? 1:p.getLocation()-(window/2);
	                	int rend = p.getLocation()+(window/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(window/2)-1;
	                	Region r = new Region(p.getGenome(), p.getChrom(), (rstart<reg.getStart()?rstart:reg.getStart()), (rend>reg.getEnd()?rend:reg.getEnd()));
	            		EnrichedFeature peak = new EnrichedFeature(r);
	                	peak.peak=p;
	                	peak.score=pval;
	                	peak.overrep=over;
	                	pset.add(peak);
	            	}
                }else if(words.length>=1){ 			//Regions only
	            	RegionParser parser = new RegionParser(gen);
	            	Region reg = parser.execute(words[0]);
	            	Point p = new Point(reg.getGenome(), reg.getChrom(), (reg.getStart()+reg.getEnd())/2);
	            	if(reg!=null){
	            		int rstart = p.getLocation()-(window/2)<1 ? 1:p.getLocation()-(window/2);
	                	int rend = p.getLocation()+(window/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(window/2)-1;
	                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
	            		EnrichedFeature peak = new EnrichedFeature(r);
	                	peak.peak=p;
                		pset.add(peak);
                	}	            	
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(pset);
	}
	
	//load towers
	public void loadTowers(String fname){
		EnrichedFeatureFileReader reader = new EnrichedFeatureFileReader(gen, fname);
		towers = reader.getFeatures();		
	}
	
	
	//Print the distribution of reads around the peaks
	//Used for calculating the fragmentation distribution
	public void printModel(String outFile){
		model.printToFile(outFile);		
	}
	
	//Convert a PairedHit
	public ReadHit convertToReadHit(Genome g, int id, PairedHit h, boolean left) {
		if(left)
			return new ReadHit(g, id, g.getChromName(h.leftChrom), (h.leftStrand ? h.leftPos : (h.leftPos-h.leftLength+1)), (h.leftStrand ? (h.leftPos+h.leftLength-1) : h.leftPos), h.leftStrand ? '+' : '-', h.weight);
		else
			return new ReadHit(g, id, g.getChromName(h.rightChrom), (h.rightStrand ? h.rightPos : (h.rightPos-h.rightLength+1)), (h.rightStrand ? (h.rightPos+h.rightLength-1) : h.rightPos), h.rightStrand ? '+' : '-', h.weight);
	}

}
