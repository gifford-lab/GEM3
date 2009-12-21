package edu.mit.csail.cgs.ewok.verbs.probers;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.binding.RegionProber;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipBayes;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.chipchip.SQLBayes;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.BayesLocator;
import edu.mit.csail.cgs.datasets.species.Genome;

public class BayesImmediateProbeGenerator implements RegionProber<Probe> {
    
    private ChipChipBayes data;
    private Genome genome;
    private double postThresh, strThresh;

    public BayesImmediateProbeGenerator(Genome g, BayesLocator loc) {
        data = loc.createObject();
        genome = g;
        postThresh = -1.0;
        strThresh = -1.0;
    }

    public BayesImmediateProbeGenerator(Genome g, BayesLocator loc, double pt, double st) {
        data = loc.createObject();
        genome = g;
        postThresh = pt;
        strThresh = st;
    }

    public BayesImmediateProbeGenerator(Genome g, ChipChipBayes d) {
        data = d;
        genome = g;
    }
    
    public Iterator<edu.mit.csail.cgs.datasets.chipchip.Probe> execute(Region r) {
        try {
            if(strThresh > 0.0 && postThresh > 0.0) { 
                data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            } else { 
                data.window(r.getChrom(), r.getStart(), r.getEnd(), postThresh, strThresh);
            }
            
            LinkedList<edu.mit.csail.cgs.datasets.chipchip.Probe> probes = 
                new LinkedList<edu.mit.csail.cgs.datasets.chipchip.Probe>();
            
            for(int i = 0; i < data.getCount(); i++) { 
                probes.addLast(createProbe(data, i, genome, r.getChrom()));
            }
            return probes.iterator();

        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private static edu.mit.csail.cgs.datasets.chipchip.Probe createProbe(ChipChipBayes data, int index, Genome genome, String chrom) { 
        edu.mit.csail.cgs.datasets.chipchip.Probe p = 
            new edu.mit.csail.cgs.datasets.chipchip.Probe(genome, chrom, data.getPos(index));
        
        double[] total = new double[2];
        
        total[0] = data.getStrength(index);
        total[1] = data.getPosterior(index);
        
        String n = data.getName();
        p.addValue(n, total);
        
        return p;
    }
}
