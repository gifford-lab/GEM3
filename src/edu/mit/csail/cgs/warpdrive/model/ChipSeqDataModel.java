/*
 * Created on Jan 11, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.warpdrive.model;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.ArrayList;

import edu.mit.csail.cgs.datasets.chippet.WeightedRunningOverlapSum;
import edu.mit.csail.cgs.datasets.chipseq.ChipSeqHit;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.general.StrandedRegion;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.ewok.verbs.MapperIterator;
import edu.mit.csail.cgs.ewok.verbs.chipseq.*;


/**
 * Warpdrive data model for ChipSeqHits.
 * If an extension is specified, then each return hit will
 * be extended by that many bp.
 */
public class ChipSeqDataModel extends WarpModel implements RegionModel, Runnable {
    
    private Expander<Region, ChipSeqHit> baseExpander, expander;
    private ChipSeqStrandedOverlapMapper mapper;
    private WeightedRunningOverlapSum totalSum;
    private WeightedRunningOverlapSum watsonSum;
    private WeightedRunningOverlapSum crickSum;
    private int extension;
    private int shift;
    private Region region;
    private boolean newinput, reloadInput, doSums, doHits;
    private ArrayList<ChipSeqHit> results;
    private ChipSeqDataProperties props;

    public ChipSeqDataModel(int ext, Expander<Region,ChipSeqHit> ex, ChipSeqStrandedOverlapMapper m ) {
        if (ext == 0) {
            expander = ex;
        } else {
            expander = new ExtendingChipSeqExpander(ext, ex, 0);
        }
        mapper = m;
        totalSum = null;
        watsonSum = null;
        crickSum = null;
        extension = ext;
        shift=0;
        baseExpander = ex;
        newinput = false;
        reloadInput = false;
        doSums = true;
        doHits = true;
        results = new ArrayList<ChipSeqHit>();
        props = new ChipSeqDataProperties();
    }
    public ChipSeqDataProperties getProperties() {return props;}

    public WeightedRunningOverlapSum getTotalRunningOverlap() { return totalSum; }
    public WeightedRunningOverlapSum getWatsonRunningOverlap() { return watsonSum; }
    public WeightedRunningOverlapSum getCrickRunningOverlap() { return crickSum; }
    
    public double getTotalMaxOverlap() { return totalSum.getMaxOverlap(); }
    public double getWatsonMaxOverlap() { return watsonSum.getMaxOverlap(); }
    public double getCrickMaxOverlap() { return crickSum.getMaxOverlap(); }
    
    public void setExtensionAndShift(int extension, int shift) { 
        if (this.extension == extension && this.shift == shift) {
            return;
        }
    	this.extension = extension;
        this.shift = shift;
        if (extension == 0 && shift == 0) {
            expander = baseExpander;
        } else {
            expander = new ExtendingChipSeqExpander(extension, baseExpander, shift);
        }
        mapper.setShift(shift);
        mapper.setExtension(extension);
    }
    
    protected void clearValues() {
        Region r = getRegion();
        totalSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
        watsonSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
        crickSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    setExtensionAndShift(getProperties().ExtendRead,getProperties().ShiftRead);
                    if (doSums) {
                        WeightedRunningOverlapSum[] sums = mapper.execute(region);
                        watsonSum = sums[mapper.POS_SUM_INDEX];
                        crickSum = sums[mapper.NEG_SUM_INDEX];
                        totalSum = sums[mapper.TOTAL_SUM_INDEX];
                    }
                    if (doHits) {
                        Iterator<ChipSeqHit> hits = expander.execute(region);
                        results.clear();
                        while (hits.hasNext()) {
                            results.add(hits.next());
                        }
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                newinput = false;
                notifyListeners();
            }
        }
    }
    public void setDoSums(boolean b) {doSums = b;}
    public void setDoHits(boolean b) {doHits = b;}
    public Iterator<ChipSeqHit> getResults() {
        return results.iterator();
    }
    protected void doReload() {
    	synchronized(this) {
    		reloadInput = true;
			this.notifyAll();
		}
    }    
}
/**
 * Wraps an Expander from Region to ChipSeqHit such that the results
 * are the reads extended by some distance.  Uses the ChipSeqHitExtender
 */
class ExtendingChipSeqExpander implements Expander<Region,ChipSeqHit> {
    
    private int extension, shift;
    private ChipSeqHitExtender extender;
    private Expander<Region,ChipSeqHit> expander;
    
    public ExtendingChipSeqExpander(int ext, Expander<Region,ChipSeqHit> exp, int shift) { 
        extension = ext;
        this.shift =shift; 
        extender = new ChipSeqHitExtender(ext, shift);
        expander = exp;
    }

    public Iterator<ChipSeqHit> execute(Region a) {
        int ns = Math.max(a.getStart() - extension, 0);
        int ne = a.getEnd() + extension;
        return new MapperIterator(extender,
                                  expander.execute(new Region(a.getGenome(), a.getChrom(), ns, ne)));
    }     
}
class ChipSeqHitExtender implements Mapper<ChipSeqHit,ChipSeqHit> {
    int distance, shift=0;
    public ChipSeqHitExtender (int d, int s) {
        distance = d; this.shift=s;
    }
    public ChipSeqHit execute(ChipSeqHit input) {
    	if(shift>0)
    		return input.shiftExtendHit(distance, shift);
    	else
    		return input.extendHit(distance);
    }
}
