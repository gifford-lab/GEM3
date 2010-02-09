package edu.mit.csail.cgs.utils.io.parsing;

public class ParseException extends Exception {
    public ParseException(String m, Throwable e) {
        super(m,e);
    }
    public ParseException(String m) {
        super(m);
    }
    public ParseException(Throwable e) {
        super(e);
    }
}
