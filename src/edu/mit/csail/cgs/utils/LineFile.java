/*
 * Created on Apr 20, 2005
 */
package edu.mit.csail.cgs.utils;
import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class LineFile {

    private Vector<String> lines;
    
    public LineFile(File f) throws IOException {
        lines = new Vector<String>();
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        while((line = br.readLine()) != null) { 
            lines.add(line);
        }
        br.close();
    }
    
    public Collection<String> getAllLines() { return lines; }
    public String getLine(int i) { return lines.get(i); }
    public int getNumLines() { return lines.size(); }
}
