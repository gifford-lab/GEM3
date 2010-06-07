/*
 * Created on Mar 3, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.util.Map;
import java.util.Collection;
import java.util.ArrayList;
import java.util.Iterator;

import edu.mit.csail.cgs.datasets.general.NamedStrandedRegion;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.types.*;

/**
 * @author tdanford
 */
public class GeneToPromoter 
    implements Mapper<Gene,NamedStrandedRegion>, SelfDescribingVerb, DefaultConstantsParameterized {

    private int upstream, downstream;
    private ArrayList<Expander<Region, ? extends Region>> dontoverlap;
    
    /** will generate promoter regions this many bases upstream and downstream of
     *  the gene's start site.
     */
    public GeneToPromoter(int up, int down) {
		upstream = up;
		downstream = down;
        dontoverlap = null;
    }
    /** will generate promoter regions this many bases upstream and downstream of
     *  the gene's start site and that do not overlap the features returned by otherfeature.
     */
    public GeneToPromoter(int up, int down, Expander<Region, ? extends Region> otherfeature) {
		upstream = up;
		downstream = down;
        dontoverlap = new ArrayList<Expander<Region, ? extends Region>>();
        dontoverlap.add(otherfeature);
    }
    public GeneToPromoter(int up, int down, Collection<Expander<Region, ? extends Region>> otherfeatures) {
		upstream = up;
		downstream = down;
        dontoverlap = new ArrayList<Expander<Region, ? extends Region>>();
        dontoverlap.addAll(otherfeatures);
    }
    
    public GeneToPromoter() { 
        upstream = 8000;
        downstream = 2000;
        dontoverlap = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public NamedStrandedRegion execute(Gene a) {
        int start, stop;
        NamedStrandedRegion output = null;
        int chrlen = a.getGenome().getChromLength(a.getChrom());
                
        switch(a.getStrand()) { 
        case '-':
            start = a.getEnd() - downstream;
            stop = a.getEnd() + upstream;
            if (start < 0) {
                start = 0;
            }
            if (stop < 0) {
                stop = 0;
            }
            if (start >= chrlen) {
                start = chrlen - 1;
            }
            if (stop >= chrlen) {
                stop = chrlen - 1;
            }
            output = new NamedStrandedRegion(a.getGenome(), a.getChrom(), start, stop, a.getID(), a.getStrand());
            break;
        default:
        case '+':
            start = a.getStart() - upstream;
            stop = a.getStart() + downstream;
            if (start < 0) {
                start = 0;
            }
            if (stop < 0) {
                stop = 0;
            }
            if (start >= chrlen) {
                start = chrlen - 1;
            }
            if (stop >= chrlen) {
                stop = chrlen - 1;
            }
            output = new NamedStrandedRegion(a.getGenome(), a.getChrom(), start, stop, a.getID(), a.getStrand());
            break;
        }
        //        System.err.println("Gene is " + a + "  -> " + a.getStart() +","+a.getEnd()+"," + a.getStrand()); 
        //        System.err.println("Intermediate output " + output + " ->  "  + start +"," + stop);
        if (dontoverlap != null) {
            start = output.getStart();
            stop = output.getEnd();
            for (Expander<Region, ? extends Region> expander : dontoverlap) {
                Iterator<? extends Region> iter = expander.execute(output);
                while (iter.hasNext()) {
                    Region other = iter.next();
                    if (other.getStart() == a.getStart() && other.getEnd() == a.getEnd()) {
                        continue;
                    }
                    if (other.overlaps(a)) {
                        continue;
                    }


                    if (a.getStrand() == '+') {
                        start = Math.min(Math.max(start, other.getEnd()), stop);
                    } else {
                        stop = Math.max(Math.min(stop, other.getStart()), start);
                    }
                    //                     System.err.println(String.format("%s %d %d -> %d %d",
                    //                                                      other.toString(), other.getStart(),other.getEnd(), start,stop));
                }
            }
            if (start != output.getStart() || stop != output.getEnd()) {
                output = new NamedStrandedRegion(output.getGenome(), output.getChrom(), start, stop, output.getName(), output.getStrand());
            }           
        }
        return output;
    }
    
    /*
     * Methods for the implementation of SelfDescribingVerb
     */

    public EchoType getInputClass() {
        return new ClassType(Gene.class);
    }

    public EchoType getOutputClass() {
        return new ClassType(NamedStrandedRegion.class);
    }
    
    private static final EchoType[] pclasses = { 
    	new ClassType(Integer.class), new ClassType(Integer.class) };
    private static final String[] pnames = { "Upstream", "Downstream" };
    
    private static final EchoType[] inputClasses = { new ClassType(Gene.class) };
    private static final String[] inputNames = { "Genes" };

    public EchoType[] getParameterClasses() {
        return pclasses;
    }

    public String[] getParameterNames() {
        return pnames;
    }
    
    public EchoType[] getInputClasses() { return inputClasses; }
    public String[] getInputNames() { return inputNames; }

    public void init(Map<String, Object> params) {
        upstream = (Integer)params.get(pnames[0]);
        downstream = (Integer)params.get(pnames[1]);
    }

    public String[] defaultConstantNames() {
        return pnames;
    }
    
    private SelfDescribingConstant[] defConsts = { new ValueWrapper(8000), new ValueWrapper(2000) };

    public SelfDescribingConstant[] defaultConstants() {
        return defConsts;
    }

}
