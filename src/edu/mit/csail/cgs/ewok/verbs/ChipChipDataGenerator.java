package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;
import edu.mit.csail.cgs.ewok.*;
import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.*;
import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.general.Region;

public class ChipChipDataGenerator implements Expander<Region,BindingEvent>{
    private ChipChipData data;
    public ChipChipDataGenerator (ChipChipData data) {
        this.data = data;
    }
    public Iterator <BindingEvent> execute(Region r) {
        ArrayList results = new ArrayList<BindingEvent>();
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            double avgs[] = new double[data.getCount()];
            for (int i = 0; i < data.getCount(); i++) {
                int numrepl = data.getReplicates(i);
                double avg = 0; int count = 0, j = 0;
                for (j = 0; j < numrepl; j++) {
                    avg += data.getRatio(i,j);
                }
                avg = avg / numrepl;
                avgs[i] = avg;
            }
            for (int i = 1; i < data.getCount() - 1; i++) {
                results.add(new BindingEvent(r.getGenome(), r.getChrom(),
                                             data.getPos(i), data.getPos(i),
                                             avgs[i],1,"Raw observation"));
            }
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();
    }
}
