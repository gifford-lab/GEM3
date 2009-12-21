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

public class MSPLazyProbeGenerator implements Expander<Region,MSPProbe>{
    
    private MSPLocator locator;
    private Genome genome;

    public MSPLazyProbeGenerator(Genome g, MSPLocator loc) {
        locator = loc;
        genome = g;
    }

    public Iterator<MSPProbe> execute(Region r) {
        try {
            ChipChipMSP data = locator.createObject();
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            //System.out.println("++++ ChipChipMSP Opened.");
            
            return new LazyIterator(data, genome, r.getChrom());
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
    }
    
    private MSPProbe createProbe(ChipChipMSP data, int index, Genome genome, String chrom) {
        
        String n = locator.name + "," + locator.version;
        MSPProbe p = 
            new MSPProbe(genome, chrom, data.getPos(index), n, data, index);
        
        return p;
    }
    
    private class LazyIterator implements Iterator<MSPProbe>, Closeable {
        
        private ChipChipMSP data;
        private Genome genome;
        private String chrom;
        private int index;
        
        private LazyIterator(ChipChipMSP d, Genome g, String chrom) { 
            data = d;
            index = 0;
            genome = g;
            this.chrom = chrom;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            if(!isClosed() && index >= data.getCount()) { 
                close();
            }
            return data != null;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public MSPProbe next() {
            if(isClosed()) { throw new IllegalStateException(); }
            if(!hasNext()) { throw new IllegalStateException(); }
            
            MSPProbe p = createProbe(data, index, genome, chrom);
            
            //System.out.println(index + ": " + p);
            
            index += 1;
            if(index >= data.getCount()) { 
                close();
            }
            
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
            if(!isClosed()) { 
                //System.out.println("---- ChipChipMSP Closed.");
                ((SQLMSP)data).close();
                data = null;
            }
        }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
         */
        public boolean isClosed() {
            return data == null;
        } 
        
    }

}
