/*
 * Author: tdanford
 * Date: Apr 30, 2008
 */
package edu.mit.csail.cgs.utils.io.parsing.alignment;

public class ElandParsingException extends Exception {
    
    private int lineNum;
    private String line;

    public ElandParsingException(int lineNum, String line) { 
        super(String.format("Error in parsing line# %d : \"%s\"", lineNum, line));
        this.lineNum = lineNum;
        this.line = line;
    }
    
    public String getLine() { return line; }
    public int getLineNum() { return lineNum; }
    
}
