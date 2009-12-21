package edu.mit.csail.cgs.ewok.verbs.probers;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.binding.RegionProber;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;

public class ChipChipImmediateProbeGenerator implements Expander<Region,Probe>, RegionProber<Probe>, Closeable {
    
    private ChipChipLocator locator;
    private ChipChipData data;
    private Genome genome;
    private String baseKey;

    public ChipChipImmediateProbeGenerator(Genome g, ChipChipLocator loc) {
        locator = loc;
        baseKey = locator.name + "," + locator.version;
        data = loc.createObject();
        genome = g;
    }

    public Iterator<edu.mit.csail.cgs.datasets.chipchip.Probe> execute(Region r) {
        if(isClosed()) { throw new IllegalStateException(); }
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            
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
    
    private edu.mit.csail.cgs.datasets.chipchip.Probe createProbe(ChipChipData data, int index, 
            Genome genome, String chrom) {
        
        edu.mit.csail.cgs.datasets.chipchip.Probe p = 
            new edu.mit.csail.cgs.datasets.chipchip.Probe(genome, chrom, data.getPos(index));
        
        double mean = 0.0;
        int reps = data.getReplicates(index);
        double[] repArray = new double[reps];
        
        double[] ipArray = new double[reps];
        double[] wceArray = new double[reps];
        
        for(int i = 0; i < data.getReplicates(index); i++) {
            double ratio = data.getRatio(index, i);
            double ip = data.getIP(index, i);
            double wce = data.getWCE(index, i);
            mean += ratio;
            repArray[i] = ratio;
            ipArray[i] = ip;
            wceArray[i] = wce;
        }
        
        if(data.getReplicates(index) > 0) { 
            mean /= (double)data.getReplicates(index);
        }
        
        String n = locator.name + "," + locator.version;

        p.addValue(getMeanKey(), mean);
        p.addValue(getRepsKey(), repArray);
        p.addValue(getWCEKey(), wceArray);
        p.addValue(getIPKey(), ipArray);
        
        return p;
    }

    public String getMeanKey() {return baseKey + ":mean";}
    public String getRepsKey() {return baseKey + ":reps";}
    public String getWCEKey() {return baseKey + ":wce";}
    public String getIPKey() {return baseKey + ":ip";}

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        ((SQLData)data).close();
        data = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return data == null;
    }
}
