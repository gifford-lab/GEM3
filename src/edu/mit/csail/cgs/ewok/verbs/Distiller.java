package edu.mit.csail.cgs.ewok.verbs;

import java.util.Iterator;

/* A Distiller maps objects of type A to objects of type B.
   Unlike a Mapper or an Expander, the Distiller may consume
   multiple objects of type A through it's execute() method before
   producing a new output of type B (available both through execute()
   when it is first produced and later through getCurrent()) */

public interface Distiller<A,B> {
    public B execute(A a);
    public B getCurrent();
    /* clears the distiller's state */
    public void reset();
}
