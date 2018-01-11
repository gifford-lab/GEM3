package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.ewok.verbs.Expander;

/** An ExpanderModel wraps an ewok Expander (something that maps an
 * input of type A to an iterator over type B).  The expander's execute method
 *  is called in response to new inputs and the results are cached
 * such that the getResults() method can be called multiple times.
 *
 * This caching is critical to the visualizer- we don't want to rerun the database query
 * every time the screen is redrawn (eg when the window is resized or uncovered) but the
 * input hasn't changed 
 */

public class ExpanderModel<IN,OUT> extends WarpModel implements Runnable {

    private boolean newinput;
    public boolean reloadInput = false;
    private Expander<IN,OUT> expander;
    private IN input;
    private ArrayList<OUT> result;

    public ExpanderModel(Expander<IN,OUT> ex) {
        super();
        expander = ex;
        newinput = false;
    }

    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!(newinput || reloadInput)) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput || reloadInput) {
            	reloadInput = false;
                try {
                    Iterator<OUT>iter = expander.execute(input);
                    result = new ArrayList<OUT>();
                    clearValues();
                    
                    while (iter.hasNext()) {
                        OUT nextValue = iter.next();
                        if(registerValue(nextValue)) { 
                            result.add(nextValue);
                        }
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                newinput = false;
                notifyListeners();
            }
        }
        
//        System.err.println("ExpanderModel run() is finishing.");
    }

    protected void setExpander(Expander<IN,OUT> expander) {
        synchronized(this) {
            if (this.expander != null && this.expander instanceof Closeable) {
                ((Closeable)this.expander).close();
            }
            this.expander = expander;
        }
    }

    public void doReload() {
    	synchronized(this) {
    		reloadInput = true;
			this.notifyAll();
		}
    }
    
    protected void clearValues() { 
    }
    
    protected boolean registerValue(OUT value) {
        return true;
    }
    
    public synchronized void setInput(IN i) {
        if (newinput == false) {
            if (!i.equals(input)) {
                input = i;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }

    public Iterator<OUT> getResults() {
        if (newinput == false && result != null) {
            return result.iterator();
        } else {
            // or should we throw some exception?
            return null;
        }
    }

    public boolean isReady() {return !newinput;}
}
