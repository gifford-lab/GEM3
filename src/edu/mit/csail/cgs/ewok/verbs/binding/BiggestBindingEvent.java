package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.binding.BindingEvent;
import edu.mit.csail.cgs.ewok.verbs.Distiller;

public class BiggestBindingEvent<X extends BindingEvent> implements Distiller<X,X> {

    private X biggest;
    public BiggestBindingEvent() {
        biggest = null;
    }
    public X execute(X input) {
        if (biggest == null || 
            input.getSize() > biggest.getSize()) {
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
