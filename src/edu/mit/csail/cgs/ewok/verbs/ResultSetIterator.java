/*
 * Created on Apr 4, 2006
 */
package edu.mit.csail.cgs.ewok.verbs;

import java.sql.*;
import java.util.Iterator;

import edu.mit.csail.cgs.utils.Closeable;

/**
 * @author tdanford
 * 
 * Adapts a java.sql.ResultSet to the Iterator interface.
 */
public class ResultSetIterator implements Iterator<ResultSet>, Closeable {
    
    private ResultSet rs;
    private boolean isReady;
    private boolean hasNext;

    public ResultSetIterator(ResultSet rs) {
        this.rs = rs;
        isReady = false;
        hasNext();
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    public boolean hasNext() {
        if(!isReady) { 
            try {
                hasNext = rs.next();
                isReady = true;
                
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
        return hasNext; 
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    public ResultSet next() {
        if(!hasNext()) { throw new IllegalStateException(); }
        isReady = false;
        return rs; 
    }

    /* (non-Javadoc)
     * @see java.util.Iterator#remove()
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }
    
    public boolean isClosed() { return rs == null; }
    
    public void close() { 
        try {
            rs.close();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        
        rs = null;
    }

}
