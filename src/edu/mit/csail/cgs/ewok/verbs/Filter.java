package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

/** a filter maps an object of type A to an object of type
 * B or to null
 */
public interface Filter<A,B> {

    public B execute(A a);
    
    public static class Compose<x,y,z> implements Filter<x,z> {
        
        private Filter<x,y> first;
        private Filter<y,z> second;
        
        public Compose(Filter<x,y> f1, Filter<y,z> f2) { 
            first = f1; 
            second = f2;
        }
        
        public z execute(x a) { 
            y firstres = first.execute(a);
            if(firstres == null) { return null; }
            return second.execute(firstres);
        }
    }
    
    public static class Identity<X> implements Filter<X,X> {
        
        public Identity() {}

        /* (non-Javadoc)
         * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(java.lang.Object)
         */
        public X execute(X a) {
            return a;
        } 
        
    }
}
