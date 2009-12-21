/*
 * Created on Mar 8, 2006
 */
package edu.mit.csail.cgs.datasets.binding;

import java.util.*;

import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class ProbedBindingEvent extends edu.mit.csail.cgs.datasets.binding.BindingEvent {
    
    private LinkedList<Probe> probes;
    
    public ProbedBindingEvent(Genome g, String c, int start, int end, 
            double size, double conf, String type, Collection<Probe> prbs) {
        super(g, c, start, end, size, conf, type);
        probes = new LinkedList<Probe>(prbs);
    }
    
    public ProbedBindingEvent(Genome g, String c, int start, int end, 
            double size, double conf, String type, Probe prb) {
        super(g, c, start, end, size, conf, type);
        probes = new LinkedList<Probe>();
        probes.addLast(prb);
    }    
    
    public ProbedBindingEvent(Genome g, String c, int start, int end, 
            double size, double conf, String type) {
        super(g, c, start, end, size, conf, type);
        probes = new LinkedList<Probe>();
    }
    
    public ProbedBindingEvent(BindingEvent evt, Collection<Probe> prbs) { 
        super(evt);
        probes = new LinkedList<Probe>(prbs);
    }

    public Collection<Probe> getProbes() { return probes; }
    
    public BindingEvent combine(BindingEvent base) {
        BindingEvent evt = super.combine(base);
        LinkedList<Probe> newProbes = probes;
        if(base instanceof ProbedBindingEvent) { 
            newProbes = new LinkedList<Probe>(probes);
            newProbes.addAll(((ProbedBindingEvent)base).getProbes());
        }
        ProbedBindingEvent pbv = new ProbedBindingEvent(evt, newProbes);
        return pbv;
    }

}
