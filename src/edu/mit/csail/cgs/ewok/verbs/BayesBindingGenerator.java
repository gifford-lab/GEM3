package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.ewok.types.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.binding.BindingExtent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.locators.BayesLocator;

public class BayesBindingGenerator<X extends Region>
	implements Expander<X,BindingExtent>, Closeable, SelfDescribingVerb, DefaultConstantsParameterized {

    private ChipChipBayes data;
    private double sizethresh, probthresh, sizeprior, probprior;
    private int totalcreated = 0;
    private boolean peaks;

    public BayesBindingGenerator (ChipChipBayes data, double probthresh, double sizethresh, boolean peaks) {
        this.data = data;
        this.sizethresh = sizethresh;
        this.probthresh = probthresh;
        this.peaks = peaks;
        sizeprior = .0;
        probprior = (probthresh < .2) ? probthresh : .2;
    }

    public BayesBindingGenerator (ChipChipBayes data, double probthresh, double sizethresh, double probprior, double sizeprior) {
        this.data = data;
        this.sizethresh = sizethresh;
        this.probthresh = probthresh;
        this.peaks = true;
        this.sizeprior = sizeprior;
        this.probprior = probprior;
    }
    
    public BayesBindingGenerator() { 
        data = null;
        sizethresh = probthresh = 0.0;
        peaks = true;
        sizeprior = probprior = 0.0;
    }

    public void init(Map<String, Object> params) {
        close();
        data = ((BayesLocator)params.get("Experiment")).createObject();
        probthresh = (Double)params.get("ProbThresh");
        sizethresh = (Double)params.get("SizeThresh");
        sizeprior = 0.0;
        probprior = (probthresh < 0.2) ? probthresh : 0.2;
    }

    public void setPeaks(boolean peaks) {
        this.peaks = peaks;
    }
    
    public boolean isClosed() { return data == null; }
    
    public void close() { 
        if(data instanceof Closeable) { 
            Closeable c = (Closeable)data;
            c.close();
        }
        data = null;
    }

    public Iterator<BindingExtent> execute(X r) {
        if (peaks) {
            return executePeaks(r);
        } else {
            ArrayList results = new ArrayList<BindingExtent>();
            try {
                data.window(r.getChrom(),
                            r.getStart(),
                            r.getEnd(),
                            probthresh, sizethresh);                
                int count = data.getCount();
                for (int i = 0; i < count; i++) {
                    int bindpos = data.getPos(i);
                    results.add(new BindingExtent(r.getGenome(),
                                                 r.getChrom(),
                                                 bindpos,
                                                 bindpos,
                                                 data.getStrength(i),
                                                 data.getPosterior(i),
                                                 "Bayes",
                                                  bindpos,
                                                  bindpos));
                }
            } catch (NotFoundException ex) {
                throw new RuntimeException("Couldn't window()" +ex,ex);
            }
            return results.iterator();
        }
    }

    /* executePeaks creates binding events as follows: first it finds
       contiguous regions in which the posterior is > probprior and
       the strength is > sizeprior.  Within this region, it requires
       the maximum posterior to be greater than probthresh and the sum
       of posterior*strength > sizethresh.  The "strength" returned in
       the BindingExtent is the sum of posterior*strength across the
       region.

    */
    public Iterator<BindingExtent> executePeaks(X r) {
        ArrayList results = new ArrayList<BindingExtent>();
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            int count = data.getCount();            
            for (int i = 0; i < count; i++) {
                if ((data.getPosterior(i) < probprior) || (data.getStrength(i) < sizeprior)) {continue;}
                double sum = 0, strsum = 0, probsum = 0, maxstr = 0, maxprob = 0;
                int j = i;
                while (j < count) {
                    if ((data.getPosterior(j) < probprior) || (data.getStrength(j) < sizeprior)) {break;}
                    if (data.getPos(j) - data.getPos(i) > 300) {break;}
                    if ((data.getPosterior(j) > maxprob) &&
                        (data.getStrength(j) > maxstr)) {
                        maxprob = data.getPosterior(j);
                        maxstr = data.getStrength(j);
                    }
                    strsum += data.getStrength(j);
                    probsum += data.getPosterior(j) * data.getStrength(j);
                    sum += data.getPosterior(j) * data.getStrength(j) * data.getPos(j);
                    j++;
                }
                if (j >= count) {
                    j = count - 1;
                }
                
                maxstr = probsum;

                if ((maxprob >= probthresh) &&
                    (maxstr >= sizethresh)) {
                    results.add(new BindingExtent(r.getGenome(),
                                                 r.getChrom(),
                                                 (int)(sum / probsum),
                                                 (int)(sum / probsum),
                                                 probsum,
                                                 maxprob,
                                                 "Bayes",
                                                  data.getPos(i),
                                                  data.getPos(j)));
                }
                i = j + 1;
            }
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();
    }

	public EchoType getInputClass() {
		return new ClassType(Region.class);
	}

	public EchoType getOutputClass() {
		return new ClassType(BindingExtent.class);
	}
	
	private static final EchoType[] paramClasses = { new ClassType(BayesLocator.class), 
		new ClassType(Double.class), new ClassType(Double.class) };
	private static final String[] paramNames = { "Experiment", "ProbThresh", "SizeThresh" };

    private static final EchoType[] inputClasses = { new ClassType(Region.class) };
	private static final String[] inputNames = { "Regions" };
    
    private static final String[] defConstNames = { "ProbThresh", "SizeThresh" };
    private static final SelfDescribingConstant[] defConsts = { new ValueWrapper(0.5), new ValueWrapper(2.0) };

	public EchoType[] getParameterClasses() {
		return paramClasses;
	}

	public String[] getParameterNames() {
		return paramNames;
	}
	
	public EchoType[] getInputClasses() { 
		return inputClasses;
	}
	
	public String[] getInputNames() { 
		return inputNames; 
	}

    public String[] defaultConstantNames() {
        return defConstNames;
    }

    public SelfDescribingConstant[] defaultConstants() {
        return defConsts;
    }

}
