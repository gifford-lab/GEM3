package edu.mit.csail.cgs.warpdrive.model;

import java.io.IOException;
import java.util.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.projects.readdb.*;
import edu.mit.csail.cgs.utils.probability.NormalDistribution;
import edu.mit.csail.cgs.utils.stats.StatUtil;

/**
 * Data model for chipseq histogram.  Separate methods for retrieving
 * plus and minus strand results
 */
public class ChipSeqHistogramModel extends WarpModel implements RegionModel, Runnable {
    
    private Client client;
    private Aggregator aggregator;
    private Map<Integer,Float> resultsPlus, resultsMinus;
    private Set<ChipSeqAlignment> alignments;
    private Set<String> ids;
    private ChipSeqHistogramProperties props;
    private double[] gaussian;	//Gaussian kernel for density estimation
    private int kernelWidth = 0;

    private Region region;
    private boolean newinput;

    public ChipSeqHistogramModel (ChipSeqAlignment a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.add(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        aggregator = new Aggregator(client);
        ids = new HashSet<String>();
        ids.add(Integer.toString(a.getDBID()));
    }
    public ChipSeqHistogramModel (Collection<ChipSeqAlignment> a) throws IOException, ClientException {
        alignments = new HashSet<ChipSeqAlignment>();
        alignments.addAll(a);
        props = new ChipSeqHistogramProperties();
        region = null;
        newinput = false;
        client = new Client();
        aggregator = new Aggregator(client);
        ids = new HashSet<String>();
        for (ChipSeqAlignment align : alignments) {
            ids.add(Integer.toString(align.getDBID()));
        }
    }    
    public ChipSeqHistogramProperties getProperties() {return props;}
    
  //pre-calculate and store the Guassian kernel prob., for efficiency
    private void initGaussianKernel(int width){
    	kernelWidth = width;
		gaussian = new double[250]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, width*width);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
    }
    
    public void clearValues() {
        resultsPlus = null;
        resultsMinus = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public Map<Integer,Float> getPlus() {return resultsPlus;}
    public Map<Integer,Float> getMinus() {return resultsMinus;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    int width = props.BinWidth;
                    int extension = props.ReadExtension;
                    // for GaussianKernel, get 1bp resolution data
                    if (props.GaussianKernelWidth!=0 && region.getWidth()<=1000){ 
                    	width = 1;
                    }
                    if (props.UseWeights) {
                        resultsPlus = aggregator.getWeightHistogram(ids,
                                                                    region.getChrom() + '+',
                                                                    region.getStart(),
                                                                    region.getEnd(),
                                                                    width,
                                                                    Float.NaN,
                                                                    extension);
                        resultsMinus = aggregator.getWeightHistogram(ids,
                                                                     region.getChrom() + '-',
                                                                     region.getStart(),
                                                                     region.getEnd(),
                                                                     width,
                                                                     Float.NaN,
                                                                     extension);
                    } else {
                        Map<Integer,Integer> tmp = aggregator.getHistogram(ids,
                                                                           region.getChrom() + '+',
                                                                           region.getStart(),
                                                                           region.getEnd(),
                                                                           width,
                                                                           Float.NaN,
                                                                           extension);
                        resultsPlus = new TreeMap<Integer,Float>();
                        for (int i : tmp.keySet()) {
                            resultsPlus.put(i, (float)tmp.get(i));
                        }
                        tmp = aggregator.getHistogram(ids,
                                                      region.getChrom() + '-',
                                                      region.getStart(),
                                                      region.getEnd(),
                                                      width,
                                                      Float.NaN,
                                                      extension);
                        resultsMinus = new TreeMap<Integer,Float>();
                        for (int i : tmp.keySet()) {
                            resultsMinus.put(i, (float)tmp.get(i));
                        }
                    }
                    // Gaussian kernel density plot
                    // update the two data hashtables with probability
                    if (props.GaussianKernelWidth!=0){
//                    	//long tic = System.currentTimeMillis();
//                    	if (gaussian==null){
//                    		this.initGaussianKernel(props.GaussianKernelWidth);
//                    	}
//                    	else if (props.GaussianKernelWidth!=kernelWidth){
//                    		this.initGaussianKernel(props.GaussianKernelWidth);
//                    	}
//                    	
//                    	int readCount = convertKernelDensity(resultsPlus);
//                    	readCount += convertKernelDensity(resultsMinus);
//                    	props.setTotalReadCount(readCount);
//                    	notifyListeners();
//                    	System.out.println((System.currentTimeMillis()-tic)/10*0.01+" sec in density esitmation");
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    // assign empty output.  This is useful because Client
                    // throws an exception for non-existant chromosomes, such
                    // as those for which there were no alignment results
                    resultsPlus = new TreeMap<Integer,Float>();
                    resultsMinus = resultsPlus;
                }
                newinput = false;
                notifyListeners();
            }
        }
        System.err.println("ChipSeqHistogram Model is closing");
        client.close();
    }                     
    // convert a weight histogram to a gaussian kernel density profile
    private int convertKernelDensity(Map<Integer,Float> results){

    	int count = results.size();
    	int min= Integer.MAX_VALUE;
    	int max= Integer.MIN_VALUE;
    	for (int pos: results.keySet()){
    		if (min>pos)
    			min = pos;
    		if (max<pos)
    			max= pos;
    	}
    	// get all reads, convert to basepair-resolution density
    	double[] profile = new double[max-min+1+100];	// add 50bp padding to the ends
    	for (int pos: results.keySet()){
    		profile[pos-min+50]=(double)results.get(pos);
    	}
    	
    	
    	double[] densities = StatUtil.symmetricKernelSmoother(profile, gaussian);
//    		StatUtil.gaussianSmoother(profile, props.GaussianKernelWidth);
    	// set density values back to the positions (only at certain resolution
    	// so that we can paint efficiently
    	int step = 1;
    	if (max-min>1024)
    		step = (max-min)/1024;
    	results.clear();
    	for (int i=min-50; i<=max+50; i+=step){
    		results.put(i, (float)densities[i-min+50]);
    	}
    	return count;
    }

}