package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

/** a Mapper is a filter whose output cannot be null, eg it
 *  maps an object of type A to an object of type B
 *
 */
public interface Mapper<A,B> extends Filter<A,B> {

    public B execute(A a);
    
    public static class Compose<X,Y,Z> implements Mapper<X,Z> {
        
        private Mapper<X,Y> first;
        private Mapper<Y,Z> second;
        
        public Compose(Mapper<X,Y> f1, Mapper<Y,Z> f2) { 
            first = f1; 
            second = f2;
        }
        
        public Z execute(X a) { 
            return second.execute(first.execute(a));
        }
    }
    
    public static class Identity<X> implements Mapper<X,X> {
        
        public Identity() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Mapper#execute(java.lang.Object)
         */
        public X execute(X a) {
            return a;
        } 
        
    }
}
