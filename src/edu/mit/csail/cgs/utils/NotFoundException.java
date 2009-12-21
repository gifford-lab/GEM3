package edu.mit.csail.cgs.utils;

/**
 * NotFoundException indicates that a resouce or target of a databse query couldn't
 * be found.
 *
 * @author <a href="mailto:arolfe@mit.edu">Alex Rolfe</a>
 * @version 1.0
 */
public class NotFoundException extends Exception {

    public NotFoundException(String s) {
        super(s);
    }
    public NotFoundException(String s, Exception e) {
        super(s,e);
    }

}
