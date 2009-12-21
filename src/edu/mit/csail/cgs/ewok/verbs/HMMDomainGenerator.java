package edu.mit.csail.cgs.ewok.verbs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ResourceBundle;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import be.ac.ulg.montefiore.run.jahmm.ViterbiCalculator;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.nouns.SimpleDomain;
import edu.mit.csail.cgs.utils.NotFoundException;

public class HMMDomainGenerator implements Expander<Region,SimpleDomain> {

	public static String[] variableNames = {"Value"};
    public static final String[] stateNames = {"unbound","mid","domain"};
    public static final int UNBOUND = 0, MID=1,DOMAIN = 2;

    public static final int numVariables = 1;
    public static final int numStates = 3;
    private Hmm<ObservationReal>  hmm;
	
	private ChipChipData ccd;
    private double[][] currenttransition;
    private double[] initialProbabilities;
    
    public HMMDomainGenerator(ChipChipLocator loc) { 
        init(new PropertyHMMParameters("edu.mit.csail.cgs.tools.binding.3state_optimized_HMM"), loc);
    }
    
    public HMMDomainGenerator(String props, ChipChipLocator loc) { 
        init(new PropertyHMMParameters(props), loc);
    }
    
    public HMMDomainGenerator(HMMParameters p, ChipChipLocator loc) { 
        init(p, loc);
    }
	
	public void init(HMMParameters p, ChipChipLocator loc) { 
		ccd = loc.createObject();
        
        currenttransition = new double[numStates][];
        initialProbabilities = p.getParams("initial_probabilities");
		
		double[] means, covars;
		
		if(loc.name.startsWith("Mm H3K27me3")) { 
			if(loc.name.indexOf("ES Stage") != -1) { 
                means = p.getParams("k27esmeans");
                covars = p.getParams("k27escovars");
                currenttransition[0] = p.getParams("k27es_unbound_transition");
                currenttransition[1] = p.getParams("k27es_mid_transition");
                currenttransition[2] = p.getParams("k27es_domain_transition");
            } else if (loc.name.indexOf("ES+2d Stage, before RA") != -1) { 
				means = p.getParams("k27es2means");
				covars = p.getParams("k27es2covars");
                currenttransition[0] = p.getParams("k27es2_unbound_transition");
                currenttransition[1] = p.getParams("k27es2_mid_transition");
                currenttransition[2] = p.getParams("k27es2_domain_transition");
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) { 
				means = p.getParams("k27es2p8hmeans");
				covars = p.getParams("k27es2p8hcovars");	
                currenttransition[0] = p.getParams("k27es2p8h_unbound_transition");
                currenttransition[1] = p.getParams("k27es2p8h_mid_transition");
                currenttransition[2] = p.getParams("k27es2p8h_domain_transition");
			} else if (loc.name.indexOf("2+1 day") != -1) { 
				means = p.getParams("k27es2p1means");
				covars = p.getParams("k27es2p1covars");
                currenttransition[0] = p.getParams("k27es2p1_unbound_transition");
                currenttransition[1] = p.getParams("k27es2p1_mid_transition");
                currenttransition[2] = p.getParams("k27es2p1_domain_transition");
			} else if (loc.name.indexOf("Olig2 Stage") != -1) { 
				means = p.getParams("k27olig2means");
				covars = p.getParams("k27olig2covars");	
                currenttransition[0] = p.getParams("k27olig2_unbound_transition");
                currenttransition[1] = p.getParams("k27olig2_mid_transition");
                currenttransition[2] = p.getParams("k27olig2_domain_transition");
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = p.getParams("k27hb9means");
				covars = p.getParams("k27hb9covars");
                currenttransition[0] = p.getParams("k27hb9_unbound_transition");
                currenttransition[1] = p.getParams("k27hb9_mid_transition");
                currenttransition[2] = p.getParams("k27hb9_domain_transition");
			} else { 
				means = p.getParams("k27esmeans");
				covars = p.getParams("k27escovars");
                currenttransition[0] = p.getParams("generic_unbound_transition");
                currenttransition[1] = p.getParams("generic_domain_transition");
			}
		} else if(loc.name.startsWith("Mm H3K79me2")) { 
			if(loc.name.indexOf("mES") != -1) { 
                means = p.getParams("k79esmeans");
                covars = p.getParams("k79escovars");
                currenttransition[0] = p.getParams("k79es_unbound_transition");
                currenttransition[1] = p.getParams("k79es_domain_transition");
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) {
				means = p.getParams("k79es2p8hmeans");
                covars = p.getParams("k79es2p8hcovars");
                currenttransition[0] = p.getParams("k79es2p8h_unbound_transition");
                currenttransition[1] = p.getParams("k79es2p8h_domain_transition"); 
			} else if (loc.name.indexOf("ES+2d Stage") != -1) {
				means = p.getParams("k79es2means");
                covars = p.getParams("k79es2covars");
                currenttransition[0] = p.getParams("k79es2_unbound_transition");
                currenttransition[1] = p.getParams("k79es2_domain_transition");
			} else if (loc.name.indexOf("2+1 day") != -1) {
				means = p.getParams("k79es2p1means");
                covars = p.getParams("k79es2p1covars");
                currenttransition[0] = p.getParams("k79es2p1_unbound_transition");
                currenttransition[1] = p.getParams("k79es2p1_domain_transition");
			} else if (loc.name.indexOf("Olig2 Stage") != -1) {
				means = p.getParams("k79olig2means");
				covars = p.getParams("k79olig2covars");		
				currenttransition[0] = p.getParams("k79olig2_unbound_transition");
                currenttransition[1] = p.getParams("k79olig2_domain_transition");
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = p.getParams("k79hb9means");
				covars = p.getParams("k79hb9covars");
				currenttransition[0] = p.getParams("k79hb9_unbound_transition");
                currenttransition[1] = p.getParams("k79hb9_domain_transition");
			} else { 
				means = p.getParams("k79esmeans");
				covars = p.getParams("k79escovars");		
				currenttransition[0] = p.getParams("generic_unbound_transition");
                currenttransition[1] = p.getParams("generic_domain_transition");
			}
		} else if(loc.name.startsWith("Mm Suz12")) { 
			if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) {
				means = p.getParams("suz12es2p8hmeans");
                covars = p.getParams("suz12es2p8hcovars");
                currenttransition[0] = p.getParams("suz12es2p8h_unbound_transition");
                currenttransition[1] = p.getParams("suz12es2p8h_domain_transition"); 
			} else if (loc.name.indexOf("ES+2d Stage") != -1) {
				means = p.getParams("suz12es2means");
                covars = p.getParams("suz12es2covars");
                currenttransition[0] = p.getParams("suz12es2_unbound_transition");
                currenttransition[1] = p.getParams("suz12es2_domain_transition");
			}	else{
				means = p.getParams("suz12es2means");
                covars = p.getParams("suz12es2covars");
                currenttransition[0] = p.getParams("suz12es2_unbound_transition");
                currenttransition[1] = p.getParams("suz12es2_domain_transition");
			}
		}else {//H3K4me3 
			if(loc.name.indexOf("ES Stage") != -1) { 
				means = p.getParams("k4esmeans");
				covars = p.getParams("k4escovars");
				currenttransition[0] = p.getParams("k4es_unbound_transition");
                currenttransition[1] = p.getParams("k4es_domain_transition");
			} else if (loc.name.indexOf("ES+2d Stage, before RA") != -1) { 
				means = p.getParams("k4es2means");
				covars = p.getParams("k4es2covars");				
				currenttransition[0] = p.getParams("k4es2_unbound_transition");
                currenttransition[1] = p.getParams("k4es2_domain_transition");
			} else if (loc.name.indexOf("ES+2d Stage, 8 hours post RA") != -1) { 
				means = p.getParams("k4es2p8hmeans");
				covars = p.getParams("k4es2p8hcovars");	
				currenttransition[0] = p.getParams("k4es2p8h_unbound_transition");
                currenttransition[1] = p.getParams("k4es2p8h_domain_transition");
			} else if (loc.name.indexOf("2+1 day") != -1) { 
				means = p.getParams("k4es2p1means");
				covars = p.getParams("k4es2p1covars");		
				currenttransition[0] = p.getParams("k4es2p1_unbound_transition");
                currenttransition[1] = p.getParams("k4es2p1_domain_transition");
			} else if (loc.name.indexOf("Olig2 Stage") != -1) { 
				means = p.getParams("k4olig2means");
				covars = p.getParams("k4olig2covars");		
				currenttransition[0] = p.getParams("k4olig2_unbound_transition");
                currenttransition[1] = p.getParams("k4olig2_domain_transition");
			} else if (loc.name.indexOf("Hb9 Stage") != -1) { 
				means = p.getParams("k4hb9means");
				covars = p.getParams("k4hb9covars");
				currenttransition[0] = p.getParams("k4hb9_unbound_transition");
                currenttransition[1] = p.getParams("k4hb9_domain_transition");
			} else { 
				means = p.getParams("k4esmeans");
				covars = p.getParams("k4escovars");
				currenttransition[0] = p.getParams("generic_unbound_transition");
                currenttransition[1] = p.getParams("generic_domain_transition");
			}
		}
		
		double[][] transition = currenttransition;
		System.out.println(String.format("Means: %f\t%f", means[0], means[1]));
		System.out.println(String.format("Variance: %f\t%f", covars[0], covars[1]));
		System.out.println(String.format("Trans1: %f\t%f", currenttransition[0][0], currenttransition[0][1]));
		System.out.println(String.format("Trans2: %f\t%f", currenttransition[1][0], currenttransition[1][1]));
		
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
	
	public HMMDomainGenerator(ChipChipData ccd, double[] means, double[] covars, double[][] transition) {
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
    
}

class BindingEventIterator implements Iterator<BindingExtent> {
    
    private Iterator<SimpleDomain> doms;
    
    public BindingEventIterator(Iterator<SimpleDomain> di) { 
        doms = di;
    }

    public boolean hasNext() {
        return doms.hasNext();
    }

    public BindingExtent next() {
        return new BindingExtent(doms.next().getBindingEvent());
    }

    public void remove() {
        doms.remove();
    } 
}

interface HMMParameters {
    public double[] getParams(String key);
}

class PropertyHMMParameters implements HMMParameters { 
    private ResourceBundle bundle;
    
    public PropertyHMMParameters(String name) {
        bundle = ResourceBundle.getBundle(name);
    }

    public double[] getParams(String key) {
        String value = bundle.getString(key);
        if(value == null) { throw new IllegalArgumentException(key); }
        String[] a = value.split(",");
        double[] params = new double[a.length];
        for(int i = 0; i < a.length; i++) { 
            params[i] = Double.parseDouble(a[i]);
        }
        return params;
    }
}

