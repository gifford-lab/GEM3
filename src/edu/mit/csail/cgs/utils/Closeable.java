/*
 * Created on May 30, 2005
 */
package edu.mit.csail.cgs.utils;

/**
 * @author tdanford
 * 
 * Classes the acquire resources (databsae connections, open files, network connections, etc)
 * should implement <code>Closeable</code>.  The <code>close</code> method allows
 * client classes to ensure that resources have been released or to release those resouces
 * as soon as they're no longer needed.
 */
public interface Closeable {
    public void close(); 
    public boolean isClosed();
}
