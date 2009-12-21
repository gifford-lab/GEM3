package edu.mit.csail.cgs.utils.database;

public class DatabaseException extends RuntimeException {

    public DatabaseException(String s) {
        super(s);
    }
    public DatabaseException(String s, Exception e) {
        super(s,e);
    }

}
