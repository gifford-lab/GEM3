/*
 * Created on Mar 27, 2006
 */
package edu.mit.csail.cgs.warpdrive;

import java.util.*;
import edu.mit.csail.cgs.ewok.verbs.*;

/**
 * @author tdanford
 */
public interface ValueReceiver<X, Y> extends Startable<X> {
    public void finish();
    public void addValue(Y value);
    
    public static class MappedWrapper<Base,Input,Output> implements ValueReceiver<Base,Input> {
        
        private ValueReceiver<Base,Output> inner;
        private Mapper<Input,Output> mapper;
        
        public MappedWrapper(ValueReceiver<Base,Output> in, Mapper<Input,Output> m) { 
            inner = in;
            mapper = m;
        }
        
        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.ValueReceiver#finish()
         */
        public void finish() { inner.finish(); }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.ValueReceiver#addValue(java.lang.Object)
         */
        public void addValue(Input value) {
            inner.addValue(mapper.execute(value));
        }
        
        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.Startable#start(java.lang.Object)
         */
        public void start(Base value) {
            inner.start(value);
        }
    }
    
    public static class ExpanderWrapper<Base,Input,Output> implements ValueReceiver<Base,Input> {
        
        private ValueReceiver<Base,Output> inner;
        private Expander<Input,Output> expander;
        
        public ExpanderWrapper(ValueReceiver<Base,Output> in, Expander<Input,Output> m) { 
            inner = in;
            expander = m;
        }
        
        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.ValueReceiver#finish()
         */
        public void finish() { inner.finish(); }

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.ValueReceiver#addValue(java.lang.Object)
         */
        public void addValue(Input value) {
            Iterator<Output> itr = expander.execute(value);
            while(itr.hasNext()) { 
                inner.addValue(itr.next());
            }
        }
        
        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.hyperdriveson.models.Startable#start(java.lang.Object)
         */
        public void start(Base value) {
            inner.start(value);
        }
    }
}
