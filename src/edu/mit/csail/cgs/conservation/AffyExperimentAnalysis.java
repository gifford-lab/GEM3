/*
 * Created on Sep 6, 2006
 */
package edu.mit.csail.cgs.conservation;

import java.util.*;
import java.io.*;
import java.text.*;

import edu.mit.csail.cgs.utils.parsing.affyexpr.*;

/**
 * @author tdanford
 */
public class AffyExperimentAnalysis {
    
    private AffyExperiment expt;

    public AffyExperimentAnalysis(AffyExperiment e) { 
        expt = e;
    }
    
    public AffyExperiment getExperiment() { return expt; }
    
    public double getArrangementPValue(BindingPartition part, Vector<String> tags) { 
        if(tags.size() != 4) { throw new IllegalArgumentException(); }
        double pvalue = 99.0;
        Vector<Double> pcts = getBlockPresentPercents(part, tags);
        Vector<Integer> counts = getBlockPresentCounts(part, tags);
        int totalCount = 0;
        for(int c : counts) { totalCount += c; }
        
        double fmin = Math.min(pcts.get(0), pcts.get(1)), fmax = Math.max(pcts.get(2), pcts.get(3));
        double fsep = fmin - fmax;
        if(fsep > 0.0) { 
            PartitionPValues.Predicate p = new PartitionPValues.ConservationPredicate(fsep);
            int[] sizes = new int[4];
            for(int i = 0; i < sizes.length; i++) { sizes[i] = part.getBlock(tags.get(i)).size(); }
            PartitionPValues.Results r = PartitionPValues.countArrangements(p, sizes, totalCount);
            pvalue = r.getFraction();
        }
        return pvalue;
    }
    
    public Vector<Integer> getBlockPresentCounts(BindingPartition part, Vector<String> tags) {
        Vector<Integer> pcts = new Vector<Integer>();
        for(String tag : tags) { 
            Vector<String> calls = lookupCalls(part.getBlock(tag));
            int present = 0;
            for(String c : calls) { present += c.equals("P") ? 1 : 0; }
            pcts.add(present);
        }
        return pcts;
    }
    
    public Vector<Double> getBlockPresentPercents(BindingPartition part, Vector<String> tags) {
        Vector<Double> pcts = new Vector<Double>();
        for(String tag : tags) { 
            Vector<String> calls = lookupCalls(part.getBlock(tag));
            int present = 0;
            for(String c : calls) { present += c.equals("P") ? 1 : 0; }
            double presentFrac = calls.size() > 0 ? (double)present / (double)calls.size() : 0.0;
            pcts.add(presentFrac);
        }
        return pcts;
    }
    
    public void doAnalysis(BindingPartition part, PrintStream ps) { 
        ps.println(part.getLongTitle() + " (" + part.getShortTitle() + ")");
        DecimalFormat sci = new DecimalFormat("0.000E0");
        DecimalFormat nf = new DecimalFormat("00.0");
        
        for(String tag : part.getBlockTags()) {
            Vector<Double> values = lookupValues(part.getBlock(tag));
            Vector<String> calls = lookupCalls(part.getBlock(tag));
            double mean = getMean(values), stddev = getStdDev(values);
            int present = 0;
            for(String c : calls) { present += c.equals("P") ? 1 : 0; }
            double presentFrac = calls.size() > 0 ? (double)present / (double)calls.size() : 0.0;
            
            ps.println("\t" + tag + " " + sci.format(mean) + " " + sci.format(stddev) + 
                    " " + present + " " + nf.format(presentFrac * 100.0));
        }
        ps.println();
    }
    
    public Vector<Double> lookupValues(Set<String> ids) { 
        Vector<Double> values = new Vector<Double>();
        for(String id : ids) { values.addAll(expt.lookupLocusIDValues(id)); }
        return values;
    }

    public Vector<String> lookupCalls(Set<String> ids) { 
        Vector<String> values = new Vector<String>();
        for(String id : ids) { values.addAll(expt.lookupLocusIDCalls(id)); }
        return values;
    }
    
    private static double getMean(Collection<Double> vals) { 
        double sum = 0.0;
        for(double v : vals) { sum += v; }
        if(vals.size() > 0) { sum /= (double)vals.size(); }
        return sum;
    }
    
    private static double getStdDev(Collection<Double> vals) { 
        double mean = getMean(vals);
        double sum = 0.0;
        for(double v : vals) {
            double dv = v - mean;
            sum += (dv * dv);
        }
        
        if(vals.size() > 0) {
            sum /= (double)vals.size(); 
            sum = Math.sqrt(sum);
        }
        
        return sum;
    }
}
