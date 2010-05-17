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
import edu.mit.csail.cgs.ewok.types.*;

/**
 * @author tdanford
 */
public class GeneToPromoter 
    implements Mapper<Gene,NamedStrandedRegion>, SelfDescribingVerb, DefaultConstantsParameterized {

    private int upstream, downstream;
    private ArrayList<RefGeneGenerator> refgene;
    
    public GeneToPromoter(int up, int down) {
		upstream = up;
		downstream = down;
        refgene = null;
    }
    public GeneToPromoter(int up, int down, RefGeneGenerator rg) {
		upstream = up;
		downstream = down;
        refgene = new ArrayList<RefGeneGenerator>();
        refgene.add(rg);
    }
    public GeneToPromoter(int up, int down, Collection<RefGeneGenerator> rg) {
		upstream = up;
		downstream = down;
        refgene = new ArrayList<RefGeneGenerator>();
        refgene.addAll(rg);
    }
    
    public GeneToPromoter() { 
        upstream = 8000;
        downstream = 2000;
        refgene = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public NamedStrandedRegion execute(Gene a) {
        int start, stop;
        NamedStrandedRegion output = null;
        switch(a.getStrand()) { 
        case '-':
            start = a.getEnd() - downstream;
            stop = a.getEnd() + upstream;
            output = new NamedStrandedRegion(a.getGenome(), a.getChrom(), start, stop, a.getID(), a.getStrand());
            break;
        default:
        case '+':
            start = a.getStart() - upstream;
            stop = a.getStart() + downstream;
            output = new NamedStrandedRegion(a.getGenome(), a.getChrom(), start, stop, a.getID(), a.getStrand());
            break;
        }
        //        System.err.println("Gene is " + a + "  -> " + a.getStart() +","+a.getEnd()+"," + a.getStrand()); 
        //        System.err.println("Intermediate output " + output + " ->  "  + start +"," + stop);
        if (refgene != null) {
            start = output.getStart();
            stop = output.getEnd();
            for (RefGeneGenerator generator : refgene) {
                Iterator<Gene> iter = generator.execute(output);
                while (iter.hasNext()) {
                    Gene other = iter.next();
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
