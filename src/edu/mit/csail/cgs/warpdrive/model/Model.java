/*
 * Created on Mar 25, 2006
 */
package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;
import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public interface Model extends EventSource<EventObject> {

    /* returns true if the Model is ready to be accessed.  This will
       be false if the model is retrieving values or otherwise modifying its internal 
       state */
    public boolean isReady();
    /* returns true iff the Model is still running.  When a model is done (no longer used),
       this should return false to allow any threads or other resources using the
       Model to clean themselves up when appropriate */
    public boolean keepRunning();
    /* tells the model that it should stop */
    public void stopRunning();
    public ModelProperties getProperties();
    public void notifyListeners();
}
