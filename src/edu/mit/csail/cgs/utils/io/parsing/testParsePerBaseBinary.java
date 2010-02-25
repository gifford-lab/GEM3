package edu.mit.csail.cgs.utils.io.parsing;

import java.io.*;
import java.util.*;


public class testParsePerBaseBinary {

    public static void main (String[] args) {
        try {
            ParsePerBaseBinary f = new ParsePerBaseBinary(args[0]);
            System.err.println("Start is " + f.getStart());
            System.err.println("Stop is " + f.getEnd());
            System.err.println("Stop is " + Integer.toHexString(f.getEnd()));
            
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

}
