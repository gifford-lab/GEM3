package edu.mit.csail.cgs.ewok.verbs.probers;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMSP;
import edu.mit.csail.cgs.datasets.chipchip.MSPProbe;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.chipchip.SQLMSP;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.*;
import edu.mit.csail.cgs.datasets.species.Genome;

public class MSPImmediateProbeGenerator<X extends Region> implements Expander<X,MSPProbe>, Closeable {
    
    private MSPLocator locator;
    private ChipChipMSP data;
    private Genome genome;

    public MSPImmediateProbeGenerator(Genome g, MSPLocator loc) {
        locator = loc;
        data = loc.createObject();
        genome = g;
        
        System.out.println("++++ ChipChipMSP Opened.");
    }

    public Iterator<MSPProbe> execute(X r) {
        if(isClosed()) { throw new IllegalStateException(); }
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            
            LinkedList<MSPProbe> probes = 
                new LinkedList<MSPProbe>();
            for(int i = 0; i < data.getCount(); i++) { 
                probes.addLast(createProbe(data, i, genome, r.getChrom()));
            }
            return probes.iterator();

        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private MSPProbe createProbe(ChipChipMSP data, int index, 
            Genome genome, String chrom) {
        
        String key = locator.name + "," + locator.version;
        MSPProbe p = 
            new MSPProbe(genome, chrom, data.getPos(index), key, data, index);
        
        return p;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        ((SQLMSP)data).close();
        data = null;
        System.out.println("---- ChipChipMSP Closed.");
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return data == null;
    }
}
