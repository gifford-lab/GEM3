package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipMLE;
import edu.mit.csail.cgs.datasets.general.Region;

public class MLEBindingGenerator implements Expander<Region,BindingEvent> {

    private ChipChipMLE mle;
    private double sizethresh, probthresh;

    public MLEBindingGenerator (ChipChipMLE mle, double probthresh, double sizethresh) {
        this.mle = mle;
        this.sizethresh = sizethresh;
        this.probthresh = probthresh;
    }

    public Iterator<BindingEvent> execute(Region r) {
        ArrayList results = new ArrayList<BindingEvent>();
        try {
            long t1 = System.currentTimeMillis();
            mle.window(r.getChrom(),
                       r.getStart(),
                       r.getEnd(),
                       sizethresh, probthresh);
            long t2 = System.currentTimeMillis();
            //            System.err.println("MLEBindingGenerator window() took " + (t2 - t1));
            int count = mle.getCount();
            for (int i = 0; i < count; i++) {
                results.add(new BindingEvent(r.getGenome(),
                                             r.getChrom(),
                                             mle.getPos(i),
                                             mle.getPos(i),
                                             mle.getSize(i),
                                             mle.getConf(i),
                                             "MLE"));
            }
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();
    }
    public void setPeaks(boolean peaks) {}
}
