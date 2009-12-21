package edu.mit.csail.cgs.ewok.verbs.probers;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.chipchip.SQLData;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.species.Genome;

public class ChipChipLazyProbeGenerator implements Expander<Region,Probe> {
    
    private ChipChipLocator locator;
    private Genome genome;

    public ChipChipLazyProbeGenerator(Genome g, ChipChipLocator loc) {
        locator = loc;
        genome = g;
    }

    public Iterator<edu.mit.csail.cgs.datasets.chipchip.Probe> execute(Region r) {
        try {
            ChipChipData data = locator.createObject();
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            return new LazyIterator(data, genome, r.getChrom());
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private edu.mit.csail.cgs.datasets.chipchip.Probe createProbe(ChipChipData data, int index, Genome genome, String chrom) { 
        edu.mit.csail.cgs.datasets.chipchip.Probe p = 
            new edu.mit.csail.cgs.datasets.chipchip.Probe(genome, chrom, data.getPos(index));
        
        double total = 0.0;
        for(int i = 0; i < data.getReplicates(index); i++) { 
            total += data.getRatio(index, i);
        }
        if(data.getReplicates(index) > 0) { 
            total /= (double)data.getReplicates(index);
        }
        
        //String n = data.getName();
        String n = locator.name + "," + locator.version;
        p.addValue(n, total);
        
        return p;
    }
    
    private class LazyIterator implements Iterator<Probe>, Closeable {
        
        private ChipChipData data;
        private Genome genome;
        private String chrom;
        private int index;
        
        private LazyIterator(ChipChipData d, Genome g, String chrom) { 
            data = d;
            index = 0;
            genome = g;
            this.chrom = chrom;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            if(data != null && index >= data.getCount()) { 
                ((SQLData)data).close();
                data = null;                
            }
            
            return data != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public edu.mit.csail.cgs.datasets.chipchip.Probe next() {
            if(isClosed()) { throw new IllegalStateException(); }
            if(!hasNext()) { throw new IllegalStateException(); }
            
            edu.mit.csail.cgs.datasets.chipchip.Probe p = createProbe(data, index, genome, chrom);
            
            //System.out.println(index + ": " + p);
            
            index += 1;
            if(index >= data.getCount()) { close(); }
            
            return p;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        }

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
}
