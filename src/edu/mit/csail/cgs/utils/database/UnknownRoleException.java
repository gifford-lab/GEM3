package edu.mit.csail.cgs.utils.database;

public class UnknownRoleException extends RuntimeException {
    public UnknownRoleException(String s) {
        super(s);
    }
    public UnknownRoleException(String s, Exception e) {
        super(s,e);
    }

}
