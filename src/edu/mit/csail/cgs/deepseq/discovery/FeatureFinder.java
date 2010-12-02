package edu.mit.csail.cgs.deepseq.discovery;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.deepseq.features.Feature;
import edu.mit.csail.cgs.deepseq.utilities.AnnotationLoader;
import edu.mit.csail.cgs.tools.utils.Args;
import edu.mit.csail.cgs.utils.ArgParser;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;

/**
 * Requires that all FeatureFinders contain an execute method that returns a list of Features. 
 * Also contains some variables common to all tools. 
 * @author shaunmahony
 *
 */
public abstract class FeatureFinder {
	protected List<Feature> signalFeatures=new ArrayList<Feature>(); //Signal channel features
	protected List<Feature> controlFeatures=new ArrayList<Feature>(); //Control channel features
	protected boolean stranded=false;
	protected Genome gen=null;
	protected double genomeLen=0;
	protected double mappableGenome=0;
	protected int seqwin = -1; //-1: all of the peak region sequence is printed, positive integer is window size otherwise
	protected String outName="out", outClosePeakName = "closePeak.out";
	protected boolean scanGenesOnly=false;
	protected boolean annotOverlapOnly=false;
	protected int maxAnnotDistance=50000;
	protected ArrayList<AnnotationLoader> geneAnnotations = new ArrayList<AnnotationLoader>();
	protected ArrayList<AnnotationLoader> otherAnnotations = new ArrayList<AnnotationLoader>();
	protected boolean dbconnected=false;
	protected double readLength=32;
    protected int maxThreads = 8;
	
	//Constructors
	public FeatureFinder(Genome g){
		gen=g;
		genomeLen = gen.getGenomeLength();
	}
	public FeatureFinder(String[] args){
		try {
			if(args.length==0){
				printError();System.exit(1);
			}
			ArgParser ap = new ArgParser(args);
			//Load genome
			if(ap.hasKey("species")){
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				if(pair != null){
					gen = pair.cdr();
					dbconnected=true;
					genomeLen = gen.getGenomeLength();
				}
			}else{
				//Make fake genome... chr lengths provided???
				if(ap.hasKey("geninfo") || ap.hasKey("g")){
					String fName = ap.hasKey("geninfo") ? ap.getKeyValue("geninfo") : ap.getKeyValue("g");
					gen = new Genome("Genome", new File(fName));
					genomeLen = gen.getGenomeLength();
				}else{
				    gen = null;
				}
			}
			mappableGenome = Args.parseDouble(args, "mappable", 0.8);
			readLength = Args.parseDouble(args,"readlen",readLength);
            maxThreads = Args.parseInteger(args,"threads",maxThreads);			
			setOutName(Args.parseString(args,"out",outName));
			setSeqwin(Args.parseInteger(args,"seqwin",seqwin));
			//Load annotations
			setAnnotOverlapOnly(Args.parseFlags(args).contains("annotoverlap"));
			setMaxAnnotDistance(Args.parseInteger(args,"maxannotdist",maxAnnotDistance));
			//Gene Annotations
			Collection<String> tfiles = Args.parseStrings(args,"transcripts");
			Collection<String> dbgenes = Args.parseStrings(args,"dbgenes");
			//Special case default
			if(dbgenes.size()==0 && tfiles.size()==0 && dbconnected){
			    String geneSource = (gen.getSpecies().equals("Saccharomyces cerevisiae") ||
						 gen.getSpecies().equals("Mycobacterium tuberculosis") )? 
				"sgdGene":"refGene";
			    geneAnnotations.add(new AnnotationLoader(gen, geneSource,geneSource, maxAnnotDistance, annotOverlapOnly));
			}
	        for(String s:dbgenes)
	        	geneAnnotations.add(new AnnotationLoader(gen, s, "refGene", maxAnnotDistance, annotOverlapOnly));
	        for(String s:tfiles)
	        	geneAnnotations.add(new AnnotationLoader(gen, s, "file", maxAnnotDistance, annotOverlapOnly));
	            //Other Annotations
	        Collection<String> namedRegions = Args.parseStrings(args,"namedregions");
	        for(String s:namedRegions)
	        	otherAnnotations.add(new AnnotationLoader(gen, s, "namedRegions", maxAnnotDistance, annotOverlapOnly));
	        Collection<String> namedStrandedRegions = Args.parseStrings(args,"namedstrandedregions");
	        for(String s:namedStrandedRegions)
	        	otherAnnotations.add(new AnnotationLoader(gen, s, "namedStrandedRegions", maxAnnotDistance, annotOverlapOnly));
	        Collection<String> namedTypedRegions = Args.parseStrings(args,"namedtypedregions");
	        for(String s:namedTypedRegions)
	        	otherAnnotations.add(new AnnotationLoader(gen, s, "namedTypedRegions", maxAnnotDistance, annotOverlapOnly));
	        Collection<String> otherAnnots = Args.parseStrings(args,"annots");
	        for(String s:otherAnnots)
	        	otherAnnotations.add(new AnnotationLoader(gen, s, "file", maxAnnotDistance, annotOverlapOnly));
	        if( Args.parseFlags(args).contains("repeatmasker"))
	        	otherAnnotations.add(new AnnotationLoader(gen, "repeatmasker", "repeatmasker", maxAnnotDistance, annotOverlapOnly));
	        
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//More options
		setStrandedFinding(Args.parseFlags(args).contains("stranded"));
        scanOnlyGenes(Args.parseFlags(args).contains("scangenesonly"));
	}
	
	
	//Accessors
	public String getOutName(){return outName;}
	
	public void setAnnotOverlapOnly(boolean a){annotOverlapOnly=a;}
	public void setMaxAnnotDistance(int x){maxAnnotDistance=x;}
	public void setStrandedFinding(boolean s){stranded=s;}
	public void scanOnlyGenes(boolean s){scanGenesOnly=s;}
	public void setSeqwin(int s){seqwin=s;}
	public void setOutName(String f){outName=f;}
	
	
	//Print the peaks to a file
	public void printFeatures(){printFeatures(true);}
	public void printFeatures(boolean signal){
		try{
			List<Feature> features = signal ? signalFeatures : controlFeatures;
			String fname = signal ? new String(outName+"_signal.peaks.txt") : new String(outName+"_control.peaks.txt"); 
			FileWriter fw = new FileWriter(fname);
			if (!features.isEmpty()){
				fw.write(features.get(0).headString());
			}
			for(Feature f : features)
				fw.write(f.toString());
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print coords for WarpDrive (Yuchun's method)
	public void printWarpCoords(){printWarpCoords(true);}
	public void printWarpCoords(boolean signal){
		try{
			List<Feature> features = signal ? signalFeatures : controlFeatures;
			FileWriter fw;
			// print the 10bp region for warpdrive file track
			if(signal){
				fw = new FileWriter(outName+"_signal_10bp.coord");
			}else{
				fw = new FileWriter(outName+"_control_10bp.coord");
			}
			for(Feature f : features) {
			  if (f.getPeak() != null) {
			    fw.write(f.getPeak().expand(5).toString()+"\n");
			  }
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the peaks to a file
	public void printGFF(){printGFF(true);}
	public void printGFF(boolean signal){
		try{
			List<Feature> features = signal ? signalFeatures : controlFeatures;
			String fname = signal ? new String(outName+"_signal.gff") : new String(outName+"_control.gff"); 
			FileWriter fw = new FileWriter(fname);
			if (!features.isEmpty()){
				fw.write("##gff-version 3\n");
			}
			for(Feature f : features)
				fw.write(f.toGFF());
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	//Print the peaks to stdout
	public void printFeatures2Stdout(){printFeatures2Stdout(true);}
	public void printFeatures2Stdout(boolean signal){
		List<Feature> features = signal ? signalFeatures : controlFeatures;
		for(Feature f : features){
			System.out.print(f.toString());
		}
	}
	public void printFeaturesWithinRange(int range){printFeaturesWithinRange(range, true);}
	public void printFeaturesWithinRange(int range, boolean signal)
	{
		List<Feature> features = signal ? signalFeatures : controlFeatures;
		Feature[] featuresArr = features.toArray(new Feature[0]);
		try{
			FileWriter fw = new FileWriter(outClosePeakName);
			boolean first=true;
			for(int i = 0; i < featuresArr.length-1; i++)
			{	
				for(int j = i+1; j < featuresArr.length; j++)
				{
					if(featuresArr[i].coords.getChrom().equals(featuresArr[j].coords.getChrom()))
					{
						if(first){
							fw.write(featuresArr[i].headString()); first=false;
						}
						if( Math.abs(featuresArr[i].getPeak().getLocation() - featuresArr[j].getPeak().getLocation()) <= range )
							fw.write(featuresArr[i].toString() + " --- " + featuresArr[j].toString());
					}
				}
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//end of printFeaturesWithinRange method
	
	public void printFeaturesWithinRange2Stdout(int range){printFeaturesWithinRange2Stdout(range, true);}
	public void printFeaturesWithinRange2Stdout(int range, boolean signal)
	{
		List<Feature> features = signal ? signalFeatures : controlFeatures;
		Feature[] featuresArr = features.toArray(new Feature[0]);
		boolean first=true;
		for(int i = 0; i < featuresArr.length-1; i++){
			for(int j = i+1; j < featuresArr.length; j++)
			{
				if(featuresArr[i].coords.getChrom().equals(featuresArr[j].coords.getChrom()))
				{
					if(first){
						System.out.println(featuresArr[i].headString()); first=false;
					}
					if( Math.abs(featuresArr[i].getPeak().getLocation() - featuresArr[j].getPeak().getLocation()) <= range )
						System.out.println(featuresArr[i].toString() + " --- " + featuresArr[j].toString());
				}
			}
		}
		
	}//end of printFeaturesWithinRange2Stdout method

	
	/* print out sequences around the peak */
	public void printPeakSequences(){printPeakSequences(true);}
	public void printPeakSequences(boolean signal){
		if(dbconnected){
			List<Feature> features = signal ? signalFeatures : controlFeatures;
			String fname = signal ? new String(outName+"_signal.seq") : new String(outName+"_control.seq");
			try {
				FileWriter fout = new FileWriter(fname);
				for(Feature f : features){
					fout.write(f.toSequence(seqwin));
				}fout.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}else{
			System.err.println("Not connected to the database, so cannot write sequences");
		}
	}
	
	/** add annotations to the enriched features. An annotation is a subclass of Region
     *  that is near or associated with the peak
     */
    protected void addRegionAnnotations(List<Feature> enriched) {
    	for(AnnotationLoader loader : otherAnnotations){
    		for (Feature peak : enriched) {
	            Region query;
	            if (annotOverlapOnly) {
	                query = peak.coords;
	            } else {
	                query = peak.coords.expand(maxAnnotDistance, maxAnnotDistance);
	            }
	            for(Region r : loader.getAnnotations(query)){
	                peak.addAnnotation(r);
	            }
	        }
    	}
    }
	//Find the closest genes to the enriched features
	protected void addClosestGenes(List<Feature> enriched){
		for(AnnotationLoader loader : geneAnnotations){
			for (Feature peak : enriched) {
                if (annotOverlapOnly) {
                    for(Gene gene : loader.getGenes(peak.coords)){
                        int overlap = gene.getOverlapSize(peak.coords);
                        if (peak.nearestGene == null || overlap > peak.distToGene) {
                            peak.nearestGene = gene;
                            peak.distToGene = overlap;
                        }
                    }
                    if (peak.nearestGene != null) {
                        peak.distToGene = peak.nearestGene.distance(peak.coords);
                    }
                } else {
                	peak.distToGene = maxAnnotDistance;
                	Region query = peak.coords.expand(maxAnnotDistance, maxAnnotDistance);
                	for(Gene gene : loader.getGenes(query)){
                        int distance = peak.getPeak().getLocation() - gene.getFivePrime();
                        if (gene.getStrand()=='-')
                        	distance = -distance;
                        if (Math.abs(distance) < Math.abs(peak.distToGene)) {
                            peak.nearestGene = gene;
                            peak.distToGene = distance;                            
                        }                            
                    }
                }
            }
		}
	}
	
	public abstract List<Feature> execute();	
	public abstract void printError();
}
