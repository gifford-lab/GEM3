/*
 * Created on Mar 20, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

/**
 * @author tdanford
 */
public class CastingMapper<A,B> implements Mapper<A,B> { 

    public CastingMapper() {
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.verbs.Filter#execute(null)
     */
    public B execute(A a) {
        return (B)a;
    }

}
