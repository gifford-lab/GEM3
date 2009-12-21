package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.ewok.verbs.Distiller;

public class MostLikelyBindingEvent<X extends BindingEvent> implements Distiller<X,X> {

    private X biggest;
    public MostLikelyBindingEvent() {
        biggest = null;
    }
    public X execute(X input) {
        if (biggest == null || 
            input.getConf() > biggest.getConf()) {
            biggest = input;
        }
        return biggest;
    }
    public X getCurrent() {
        return biggest;
    }
    public void reset() {
        biggest = null;
    }
}
