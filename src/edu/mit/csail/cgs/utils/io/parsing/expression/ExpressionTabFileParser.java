/*
 * Created on Sep 19, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.expression;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class ExpressionTabFileParser implements Iterator<ExprProbeValue> {
    
    private String exptName, speciesName;
    private BufferedReader br;
    private ExprProbeValue nextValue;
    private int index;
    
    public ExpressionTabFileParser(File f) throws IOException {
        index = 1;
        br = new BufferedReader(new FileReader(f));
        String firstLine = br.readLine();
        String[] array = firstLine.split("\\t");
        int count = Integer.parseInt(array[0]);
        if(count != 1) { 
        	throw new IllegalArgumentException("Count: " + count + "?"); 
        }
        speciesName = array[1];
        exptName = br.readLine();
		parseNextLine();
    }

    public ExpressionTabFileParser(File f, int index) throws IOException {
        this.index = index;
        br = new BufferedReader(new FileReader(f));
        String firstLine = br.readLine();
        String[] array = firstLine.split("\\t");
        int count = Integer.parseInt(array[0]);
        if(count != 1) { 
        	throw new IllegalArgumentException("Count: " + count + "?"); 
        }
        speciesName = array[1];
        exptName = br.readLine();
		parseNextLine();
    }
    
    public String getExptName() { return exptName; }
    public String getSpeciesName() { return speciesName; }

    private void parseNextLine() {
        nextValue = null;
        try { 
            String line = br.readLine();
			System.out.println("Expr line: " + line);
            if(line != null) { 
                String[] array = line.split("[\\s]+");
                String probe = array[0];
                double value = Double.parseDouble(array[index]);
                nextValue = new ExprProbeValue(exptName, probe, value);
            } else { 
                br.close();
            }
        } catch(IOException ie) { 
            ie.printStackTrace();
        }
    }

    public boolean hasNext() {
        return nextValue != null;
    }


    public ExprProbeValue next() {
        ExprProbeValue sv = nextValue;
        parseNextLine();
        return sv;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
    public String getDescriptor() {
        return "Expression Experiment: " + exptName;
    }
    
}
