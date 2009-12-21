package edu.mit.csail.cgs.ewok.verbs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchScaledLearner;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.species.Organism;
import edu.mit.csail.cgs.ewok.nouns.SimpleDomain;
import edu.mit.csail.cgs.utils.NotFoundException;

public class HMMParameterEstimator implements Expander<Region,SimpleDomain>{
	public static String[] variableNames = {"Value"};
    public static final String[] stateNames = {"unbound","mid", "domain"};
    public static final int UNBOUND = 0, MID=1, DOMAIN = 2;
    
    
    /**Params for H3K27me3 ES*/
    public static final double[] k27esmeans = {0.5,	//unbound
        10};	//domain
    public static final double[] k27escovars = {5,
           10};
    /**Params for H3K27me3 ES2*/
    public static final double[] k27es2means = {0.5,	//unbound
        10};	//domain
    public static final double[] k27es2covars = {1,
           20};
    /**Params for H3K27me3 ES2 + 8 hours*/
    public static final double[] k27es2p8hmeans = {0.5,	//unbound
        											10};	//domain
    public static final double[] k27es2p8hcovars = {2,
    												10};
    /**Params for H3K27me3 ES2 + 1day*/
    public static final double[] k27es2p1means = {0.5,	//unbound
        10};	//domain
    public static final double[] k27es2p1covars = {1,
           10};
    /**Params for H3K27me3 Olig2*/
    public static final double[] k27olig2means = {0.5,	//unbound
                                            10};	//domain
    public static final double[] k27olig2covars = {1,
                                               10};
    /**Params for H3K27me3 Hb9*/
    public static final double[] k27hb9means = {0.5,	//unbound
        10};	//domain
    public static final double[] k27hb9covars = {1,
           10};

    /**H3K4me3 params**/
    public static final double[] k4esmeans = {1.5,	//unbound
    											20.0};	//domain
	public static final double[] k4escovars = {0.5,
												20.0};
    public static final double[] k4es2means = {1.0,	//unbound
        										10};	//domain
    public static final double[] k4es2covars = {1.0,
    											10.0};
    public static final double[] k4es2p8hmeans = {1.0,	//unbound
        											15};	//domain
    public static final double[] k4es2p8hcovars = {1.0,
    												10.0};
    public static final double[] k4es2p1means = {0.8967,	//unbound
                                            20.0};	//domain
    public static final double[] k4es2p1covars = {2.0,
                                               10.0};
    public static final double[] k4olig2means = {0.5,	//unbound
        10};	//domain
    public static final double[] k4olig2covars = {1,
           10};
    public static final double[] k4hb9means = {0.5,	//unbound
    											10.0};	//domain
    public static final double[] k4hb9covars = {0.0385,
    											4.9176};

    /**H3K79me2 params**/
    public static final double[] k79esmeans = {1.0,	//unbound
    											20.0};	//domain
    public static final double[] k79escovars = {0.5,
    											10};
    public static final double[] k79es2means = {1.0,	//unbound
												20.0};	//domain
    public static final double[] k79es2covars = {0.5,
												10};
    public static final double[] k79es2p8hmeans = {1.0,	//unbound
    											20.0};	//domain
    public static final double[] k79es2p8hcovars = {0.5,
												10};
    public static final double[] k79es2p1means = {1.0,	//unbound
												20.0};	//domain
    public static final double[] k79es2p1covars = {0.5,
    											10};
    public static final double[] k79olig2means = {1.0,	//unbound
    											20.0};	//domain
    public static final double[] k79olig2covars = {0.5,
												10};
	public static final double[] k79hb9means = {1.0,	//unbound
												20.0};	//domain
	public static final double[] k79hb9covars = {0.5,
												10};
    public static final double[][] k79estransition = {{.995, .005},  // unbnd
        												{.005, .995}};  // domain
    
    
    public static double[] initialProbabilities = {0.5,0.25, 0.25};
    public static final double[][] generictransition = {{.9, .1, .1},  // unbnd
    													{.1, .9, .1},  //mid
        												{.1, .1, .9}};  // domain
    public static final double[] genericMeans = {1.0,	//unbound
    											5.0, //mid
												20.0};	//domain
    public static final double[] genericCovars = {0.5, 2, 10};
    
    public static final int numVariables = 1;
    public static final int numStates = 3;
    private Hmm<ObservationReal>  hmm;
	
	private ChipChipData ccd;
	public static int numIter=20;
	
	public HMMParameterEstimator(ChipChipLocator loc) { 
		ccd = loc.createObject();
		
		double[] means=genericMeans, covars=genericCovars;
		
		/*if(loc.name.startsWith("Mm H3K27me3")) { 
			if(loc.name.indexOf("ES Stage") != -1) { 
				means = k27esmeans;
				covars = k27escovars;
			} else if (loc.name.indexOf("ES+2d Stage, before RA") != -1) { 
				means = k27es2means;
				covars = k27es2covars;				
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) { 
				means = k27es2p8hmeans;
				covars = k27es2p8hcovars;			
			} else if (loc.name.indexOf("2+1 day") != -1) { 
				means = k27es2p1means;
				covars = k27es2p1covars;				
			} else if (loc.name.indexOf("Olig2 Stage") != -1) { 
				means = k27olig2means;
				covars = k27olig2covars;		
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = k27hb9means;
				covars = k27hb9covars;
			} else { 
				means = k27es2means;
				covars = k27es2covars;				
			}
		} else if(loc.name.startsWith("Mm H3K4me3")){ 
			if(loc.name.indexOf("ES Stage") != -1) { 
				means = k4esmeans;
				covars = k4escovars;
			} else if (loc.name.indexOf("ES+2d Stage, before RA") != -1) { 
				means = k4es2means;
				covars = k4es2covars;
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) { 
				means = k4es2p8hmeans;
				covars = k4es2p8hcovars;			
			} else if (loc.name.indexOf("2+1 day") != -1) { 
				means = k4es2p1means;
				covars = k4es2p1covars;				
			} else if (loc.name.indexOf("Olig2 Stage") != -1) { 
				means = k4olig2means;
				covars = k4olig2covars;		
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = k4hb9means;
				covars = k4hb9covars;
			} else { 
				means = k4es2means;
				covars = k4es2covars;				
			}
		} else if(loc.name.startsWith("Mm H3K79me2")){ 
			if(loc.name.indexOf("mES") != -1) { 
				means = k79esmeans;
				covars = k79escovars;
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) {
				means = k79es2p8hmeans;
				covars = k79es2p8hcovars;		
			} else if (loc.name.indexOf("ES+2d Stage") != -1) {
				means = k79es2means;
				covars = k79es2covars;
			} else if (loc.name.indexOf("2+1 day") != -1) {
				means = k79es2p1means;
				covars = k79es2p1covars;
			}else if (loc.name.indexOf("Olig2 Stage") != -1) { 
				means = k79olig2means;
				covars = k79olig2covars;		
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = k79hb9means;
				covars = k79hb9covars;
			} else { 
				means = k79esmeans;
				covars = k79escovars;				
			}
		}else {
			//Nothing yet...
			means = k27es2means;
			covars = k27es2covars;
		}*/
		
		double[][] transition = generictransition;
		
		ArrayList<Opdf<ObservationReal>> pdfs = new ArrayList<Opdf<ObservationReal>>();
        for (int i = 0; i < numStates; i++) {
            pdfs.add(new OpdfGaussian(means[i],
                                           covars[i]));            
        }
        hmm = new Hmm<ObservationReal>(initialProbabilities,
                                         transition,
                                         pdfs);
        for (int i = 0; i < numStates; i++) {
            System.err.println("pi[" + i + "]=" + hmm.getPi(i));
        } 
		
		this.ccd = ccd;
	}
	
	public HMMParameterEstimator(ChipChipData ccd, double[] means, double[] covars, double[][] transition) {
		ArrayList<Opdf<ObservationReal>> pdfs = new ArrayList<Opdf<ObservationReal>>();
        for (int i = 0; i < 2; i++) {
            pdfs.add(new OpdfGaussian(means[i],
                                           covars[i]));            
        }
        hmm = new Hmm<ObservationReal>(initialProbabilities,
                                         transition,
                                         pdfs);
        for (int i = 0; i < numStates; i++) {
            System.err.println("pi[" + i + "]=" + hmm.getPi(i));
        } 
		
		this.ccd = ccd;
	}
	
	/*
	 * Added 12/14/2007. CYY.
	 */
	public ArrayList<ArrayList<ObservationReal>> getObservations(List<Dataset> data) {
        ArrayList<ArrayList<ObservationReal>> observations = new ArrayList<ArrayList<ObservationReal>>();        
        for(int i=0; i<data.size(); i++) {
        	ArrayList<ObservationReal> obs = new ArrayList<ObservationReal>();
	        for (int j=0; j<data.get(i).getCount(); j++) {
	        	Datapoint point = data.get(i).getPoint(j);
	        	double val = point.getValue();
	        	if(Double.isNaN(val)||Double.isInfinite(val)){
	        		val=1.0;
	        	}
	        	obs.add(new ObservationReal(val));
	        }
	        observations.add(obs);
        }
        return observations;
    }
	
	public int[] getStates(Dataset data) {
        ArrayList<ObservationReal> observations = new ArrayList<ObservationReal>();
        for (int i=0; i<data.getCount(); i++) {
        	Datapoint point = data.getPoint(i);
        	double val = point.getValue();
        	if(Double.isNaN(val)||Double.isInfinite(val)){
        		val=1.0;
        	}
        	observations.add(new ObservationReal(val));
        }
        ViterbiCalculator vc = new ViterbiCalculator(observations,hmm);
        int[] states = vc.stateSequence();

        return states;
    }
    
    public List<SimpleDomain> getDomains(Dataset data) {
    	ArrayList<SimpleDomain> regions = new ArrayList<SimpleDomain>();
    	int[] states = getStates(data);
    	int state = states[0];
    	int start = 0;
    	int i;
    	for (i=0; i<states.length; i++) {
    		if (states[i]==DOMAIN && state != DOMAIN) {
    			start = i;
    		}
    		if (states[i]!=DOMAIN && state==DOMAIN) {
    			regions.add(new SimpleDomain(data.getGenome(), data.getChrom(), data.getPoint(start).getPos(), data.getPoint(i-1).getPos()));
    		}
    		state = states[i];
    	}
    	if (state==DOMAIN) {
    		regions.add(new SimpleDomain(data.getGenome(), data.getChrom(), data.getPoint(start).getPos(), data.getPoint(i-1).getPos()));
    	}
    	return regions;
    }
    
    public Iterator<SimpleDomain> execute(Region a) {
    	//make wellTiledRegionGenerator
    	File file = new File("C://Shaun/PPG/Datafiles/well_tiled_regions.3col.txt");
    	Iterator<Region> regions = new WellTiledRegionGenerator(file);
    	//for loop to get each region from the welltiledregiongenerator
    	ArrayList<Region> contiguousAllRegions = new ArrayList<Region>();	
    	while(regions.hasNext()){
    		contiguousAllRegions.addAll(contiguousSubRegions(regions.next()));
    	}
    	ArrayList<Dataset> wholeDatasets = new ArrayList<Dataset>();
		for (Region r : contiguousAllRegions) {
        	Dataset dataset = datasetFromContiguousRegion(r);
        	if (dataset.getCount()==0 || dataset.getCount()==1) {
//        		System.err.println(r.getChrom()+" "+r.getStart()+" "+r.getEnd());
        	} else {
        		wholeDatasets.add(dataset);
        	}
        }
		//END OF FOR LOOP
		ArrayList<ArrayList<ObservationReal>> observations = new ArrayList<ArrayList<ObservationReal>>();
		observations.addAll(getObservations(wholeDatasets));
		System.out.println("num of observation sequences: "+observations.size());
//		FileOutputStream out; // declare a file output object
//        PrintStream p; // declare a print stream object
//
//        try
//        {
//                out = new FileOutputStream("D://MIT/rotations/Gifford/observations.txt");
//
//                // Connect print stream to the output stream
//                p = new PrintStream( out );
//                
//                for(int i=0; i<observations.size(); i++){
//                	p.println((i+1)+": "+observations.get(i).toString());
//                	for(int j=0; j<observations.get(i).size(); j++){
//                		if(observations.get(i).get(j).value>10 || observations.get(i).get(j).value==0){
//                			System.out.println((i+1)+":"+(j+1)+": "+observations.get(i).get(j).value);
//                		}
//                	}
//                }
//
//                p.close();
//        }
//        catch (Exception e)
//        {
//                System.err.println ("Error writing to file");
//        }
		
//		ArrayList<ArrayList<ObservationReal>> subObservations = new ArrayList<ArrayList<ObservationReal>>();
//		for (int i=159; i<160; i++){
//			subObservations.add(observations.get(i));
//		}
		System.out.println("The initial HMM parameters:");
		System.out.println(hmm.toString());
		BaumWelchScaledLearner bwsl = new BaumWelchScaledLearner();
		bwsl.setNbIterations(numIter);
		hmm = bwsl.learn(hmm, observations);		
		int numIteration =bwsl.getNbIterations();
		System.out.println("num of interation: "+numIteration);
		System.out.println("The optimized HMM parameters:");
		System.out.println(hmm.toString());
		
		ArrayList<Region> contiguousRegions = new ArrayList<Region>();
		contiguousRegions.addAll(contiguousSubRegions(a));
		ArrayList<SimpleDomain> domains = new ArrayList<SimpleDomain>();
		for (Region r : contiguousRegions) {
        	Dataset dataset = datasetFromContiguousRegion(r);
        	if (dataset.getCount()==0) {
//        		System.err.println(r.getChrom()+" "+r.getStart()+" "+r.getEnd());
        	} else {
        		domains.addAll(getDomains(dataset));
        	}
        }
		return domains.iterator();
	}
	
	public List<Region> contiguousSubRegions(Region r) {
    	try {
			ccd.window(r.getChrom(), r.getStart(), r.getEnd());
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
//		System.err.println(r.getChrom()+" count: "+ccd.getCount());
    	ArrayList<Region> regionList = new ArrayList<Region>();
    	if (ccd.getCount()==0) {
//    		System.err.println(r.getChrom()+" empty");
    		return regionList;
    	}
    	int start = 0;
    	int end = 1;
    	while (end+1 < ccd.getCount()) {
    		if ((ccd.getPos(end+1) - ccd.getPos(end)) > (ccd.getPos(end) - ccd.getPos(end-1) + 1000)) {
    			if (end-start > 1) {
    				regionList.add(new Region(r.getGenome(), r.getChrom(), ccd.getPos(start)-1, ccd.getPos(end)+1));
    			}
    			start = end + 1;
    			end = start;
    		}
    		end++;
    	}
    	if (end-start > 1) {
    		regionList.add(new Region(r.getGenome(), r.getChrom(), ccd.getPos(start)-1, ccd.getPos(end)+1));
    	}
//    	System.err.println(r.getChrom()+" contig regions "+regionList.size());
    	return regionList;
    }
	
    public Dataset datasetFromContiguousRegion(Region r) {
    	try {
			ccd.window(r.getChrom(), r.getStart(), r.getEnd());
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
    	Dataset data = new Dataset(r.getGenome(), r.getChrom(), r);
    	for (int i=0; i<ccd.getCount(); i++) {
    		double[] reps = new double[ccd.getReplicates(i)];
			for (int j=0; j<ccd.getReplicates(i); j++) {
				reps[j] = ccd.getRatio(i, j);
			}
    		data.addDatapoint(new Datapoint(ccd.getPos(i), median(reps)));
    	}
    	return data;
    }
    
    public static double median(double[] values) {
		Arrays.sort(values);
		if (values.length % 2 == 0) {
			return (values[values.length/2-1] + values[values.length/2]) / 2.0;
		} else {
			return values[(int)Math.floor(((double)values.length)/2.0)];
		}
	}
    

    private class Datapoint implements Comparable<Datapoint> {

    	private int pos;
    	private double value;
    	
    	public Datapoint(int pos, double value) {
    		this.pos = pos;
    		this.value = value;
    	}
    	
    	public int getPos() {
    		return pos;
    	}

    	public double getValue() {
    		return value;
    	}



    	public int compareTo(Datapoint arg0) {
    		return this.pos - arg0.pos;
    	}
    	
    	public boolean equals(Object o) {
    		if (o instanceof Datapoint) {
    			return this.pos == ((Datapoint)o).pos;
    		} else {
    			return false;
    		}
    	}

    }
    
    private class Dataset {

    	private ArrayList<Datapoint> datapoints;
    	private Genome g;
    	private String chrom;
    	private Region r;
    	
    	public Dataset(Genome g, String chrom, Region r) {
    		datapoints = new ArrayList<Datapoint>();
    		this.g = g;
    		this.chrom = chrom;
    		this.r = r;
    	}
    	
    	public void addDatapoint(Datapoint d) {
            datapoints.add(d);
        }
    	
    	public void sort() {Collections.sort(datapoints);}
    	
    	public int getCount() {return datapoints.size();}
        public Datapoint getPoint(int i) {return datapoints.get(i);}
        public int getPointIndex(int p) {return Collections.binarySearch(datapoints,new Datapoint(p,0));}
    	
        public Genome getGenome() {
        	return g;
        }
        
        public String getChrom() {
        	return chrom;
        }
        
        public Region getRegion() {
        	return r;
        }
    }

    
    public static class BindingEventWrapper implements Expander<Region,BindingExtent> {
        
        private HMMDomainGenerator gen;
        
        public BindingEventWrapper(HMMDomainGenerator g) { 
            gen = g;
        }

        public Iterator<BindingExtent> execute(Region a) {
            Iterator<SimpleDomain> doms = gen.execute(a);
            return new BindingEventIterator(doms);
        } 
    }
    
    public static void main(String[] args) { 
    	Genome mm8;
		try {
			mm8 = Organism.findGenome("mm8");
//			String name = "Mm H3K27me3:HBG3:ES Stage vs H3:HBG3:ES Stage";
//			String name = "Mm H3K27me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA";
//			String name = "Mm H3K27me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA";
//			String name = "Mm H3K27me3:HBG3:2+1 day vs H3:HBG3:2+1 day";
//			String name = "Mm H3K27me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage";
			String name = "Mm H3K27me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage";

//			String name = "Mm H3K4me3:HBG3:ES Stage vs H3:HBG3:ES Stage";
//			String name = "Mm H3K4me3:HBG3:ES+2d Stage, before RA vs H3:HBG3:ES+2d Stage, before RA";
//			String name = "Mm H3K4me3:HBG3:2+1 day vs H3:HBG3:2+1 day";
//			String name = "Mm H3K4me3:HBG3:ES+2d Stage, 8 hours post RA vs H3:HBG3:ES+2d Stage, 8 hours post RA";	
//			String name = "Mm H3K4me3:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage";
//			String name = "Mm H3K4me3:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage";
			
//			String name="Mm H3K79me2:HBG3:mES vs H3:HBG3:mES";
//			String name="Mm H3K79me2:HBG3:ES+2d Stage vs WCE:HBG3:ES+2d Stage";
//			String name="Mm H3K79me2:HBG3:ES+2d Stage, 8 hours post RA vs WCE:HBG3:ES+2d Stage, 8 hours post RA";
//			String name="Mm H3K79me2:HBG3:2+1 day vs WCE:HBG3:2+1 day";
//			String name="Mm H3K79me2:HBG3:Olig2 Stage vs H3:HBG3:Olig2 Stage";
//			String name="Mm H3K79me2:HBG3:Hb9 Stage vs H3:HBG3:Hb9 Stage";
			
			numIter = 50;
			
			String version = "median linefit, quantile norm";
	    	ChipChipLocator loc = new ChipChipLocator(mm8, name, version);
	    	Expander<Region,SimpleDomain> domainCaller = new HMMParameterEstimator(loc);
	    	//make a fake region as the argument of execute(Region a).
	    	File file = new File("C://Shaun/PPG/Datafiles/well_tiled_regions.3col.txt");
	    	Iterator<Region> regions = new WellTiledRegionGenerator(file);
	    	Iterator<SimpleDomain> domains = domainCaller.execute(regions.next());
//	    	System.out.println(domains.next().toString());
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
    }
}
