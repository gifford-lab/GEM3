/*
 * Created on Feb 20, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.ewok.verbs;

public class MapperCastWrapper<A,B,C extends B> implements Mapper<A,B> {
    
    private Mapper<A,C> internalMapper;
    
    public MapperCastWrapper(Mapper<A,C> i) { internalMapper = i; }

    public B execute(A a) {
        return (B)(internalMapper.execute(a));
    }
}
