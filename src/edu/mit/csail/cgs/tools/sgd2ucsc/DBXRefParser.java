/*
 * Created on Jan 29, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.tools.sgd2ucsc;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.utils.iterators.SingleIterator;

/*
 * Columns are:
 * 1) DBXREF ID 
 * 2) DBXREF ID source
 * 3) DBXREF ID type
 * 4) S. cerevisiae feature name
 * 5) SGDID
 */

/**
 * @author tdanford
 */
public class DBXRefParser {
    
    public static void main(String[] args) { 
        File input = args.length > 0 ? new File(args[0]) : 
            new File("C:\\Documents and Settings\\tdanford\\Desktop\\dbxref.tab");
        DBXRefParser parser = new DBXRefParser(input);
        
    }
    
    public LinkedList<XRefLine> lines;
    
    public DBXRefParser(File f) { 
        lines = new LinkedList<XRefLine>();
        SingleIterator<File> fitr = new SingleIterator<File>(f);
        FileLineExpander lineExp = new FileLineExpander();
        Iterator<String> litr = new ExpanderIterator<File,String>(lineExp, fitr);
        while(litr.hasNext()) { 
            String line = litr.next();
            XRefLine xref = new XRefLine(line);
            lines.addLast(xref);
        }
        
        System.out.println("DBXRef File: " + lines.size() + " entries.");
    }

    public static class XRefLine { 
        public String id, source, type, scFeature, sgdID;
        public XRefLine(String line) { 
            String[] array = line.split("\\t");
            id = array[0]; 
            source = array[1];
            type = array[2];
            scFeature = array[3];
            sgdID = array[4];
        }
    }
}
