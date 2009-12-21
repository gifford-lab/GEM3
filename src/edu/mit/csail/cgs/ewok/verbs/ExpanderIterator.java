package edu.mit.csail.cgs.ewok.verbs;

import java.util.*;

import edu.mit.csail.cgs.utils.Closeable;

/* An ExpanderIterator takes an Expander<A,B> and an Iterator<B>
   and returns an Iterator<B> that is the result of applying
   the Expander to each element of the input iterator.
 */

public class ExpanderIterator<A,B> implements Iterator<B>, Closeable {
    private Expander<A,B> mapper;
    private Iterator<A> input;
    private Iterator<B> curiter;
    private B nextval;
    private boolean first;

    public ExpanderIterator(Expander<A,B> mapper, Iterator<A> input) {
        this.mapper = mapper;
        this.input = input;
        curiter = null;
        first = true;
        nextval = null;
    }

    private void getNext() {
    	nextval = null;
    	
    	while((curiter == null || !curiter.hasNext()) && input.hasNext()) { 
    		if(curiter != null && (curiter instanceof Closeable)) { 
    			((Closeable)curiter).close();
    		}
    		curiter = mapper.execute(input.next());
    	}
    	
    	if(curiter != null && curiter.hasNext()) { 
    		nextval = curiter.next();
    	}
    	
    	/*
        if (curiter != null &&
            curiter.hasNext()) {
        	
            nextval = curiter.next();
            //            System.err.println(" nextval="+nextval);
            return;
        } else if (curiter != null && curiter instanceof Closeable) {
            ((Closeable)curiter).close();
        }

        if (input.hasNext()) {
            A a = input.next();
            curiter = mapper.execute(a);
            //            System.err.println("ExpanderIterator just got its next input: " + a + " and curiter is " + curiter);
            //            System.err.println("  input.hasNext is " + input.hasNext() + " and curiter.hasNext is " + curiter.hasNext());
            getNext();
        } else {
            nextval = null;
            //            System.err.println(" nextval=null because input.hasNext == false");
        }
        */
    }

    public boolean hasNext() {
        if (first) {getNext();}
        first = false;
        return (nextval != null);
    }

    public B next() {
        B temp = nextval;
        getNext();
        return temp;
    }

    public void remove() {
        throw new UnsupportedOperationException("Can't remove from a ExpanderIterator");
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#close()
     */
    public void close() {
        if(input instanceof Closeable) {
            Closeable c = (Closeable)input;
            if(!c.isClosed()) { c.close(); }
        }
        
        /*
        if(mapper instanceof Closeable) {
            Closeable c = (Closeable)mapper;
            if(!c.isClosed()) { c.close(); }
        }
        */
        
        if(curiter instanceof Closeable) {
            Closeable c = (Closeable)curiter;
            if(!c.isClosed()) { c.close(); }
        }

        input = null;
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
     */
    public boolean isClosed() {
        return input==null;
    }


}
