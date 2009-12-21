package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

import edu.mit.csail.cgs.utils.Closeable;

/* An Expander<A,B> maps an object of type A to a set of objects of type
   B (represented as an Iterator<B>) */

public interface Expander<A,B> {

    public Iterator<B> execute(A a);
    
    public static class Compose<X,Y,Z> implements Expander<X,Z> {
        
        private Expander<X,Y> first;
        private Expander<Y,Z> second;
        
        public Compose(Expander<X,Y> f1, Expander<Y,Z> f2) { 
            first = f1; 
            second = f2;
        }
        
        public Iterator<Z> execute(X a) {
            return new LazyComposeIterator<X,Y,Z>(first.execute(a), second);
        }
    }
    
    public static class LazyComposeIterator<X,Y,Z> implements Iterator<Z> {
        
        private Iterator<Y> base;
        private Iterator<Z> expanded;
        private Expander<Y,Z> expander;
        
        public LazyComposeIterator(Iterator<Y> itr, Expander<Y,Z> b) { 
            base = itr;
            expander = b;
            expanded = null;
            findNext();
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#hasNext()
         */
        public boolean hasNext() {
            return expanded != null && expanded.hasNext();
        }
        
        private void findNext() {             
            while(base.hasNext() && (expanded == null || !expanded.hasNext())) {   
                expanded = expander.execute(base.next());
            }
            if(!base.hasNext()) { 
                expanded = null;
            }
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#next()
         */
        public Z next() {
            Z val = expanded.next();
            findNext();
            return val;
        }

        /* (non-Javadoc)
         * @see java.util.Iterator#remove()
         */
        public void remove() {
            throw new UnsupportedOperationException();
        } 
        
    }
}
