/*
 * Created on Mar 16, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing;

import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class ParsePerBaseBinary extends ChipChipFile {
    
    public static final long sHeaderSize = 8;
    public static final long sElementSize = 8;

    private RandomAccessFile raFile;
    private int start, end;
    
    public ParsePerBaseBinary(String filename) throws IOException {
        super(filename, "bp_dist");
        openFile();
    }
    
    private void openFile() throws IOException {
        raFile = new RandomAccessFile(new File(super.getBaseFilename()), "r");
        start = Integer.reverseBytes(raFile.readInt());
        end = Integer.reverseBytes(raFile.readInt());                                   
		System.err.println("--> " + super.getBaseFilename() + ": " + start + 
		"," + end + " : " + raFile.length());
    }
    
    public void close() throws IOException { 
        raFile.close();
        raFile = null;
    }

    public int getStart() { return start; }
    public int  getEnd() { return end; }
    
    public double getValue(int offset) throws IOException { 
        if(offset < start || offset > end) { 
            throw new IllegalArgumentException();
        }
        raFile.seek(sHeaderSize + ((long)(offset - start) * sElementSize));
        return Double.longBitsToDouble(Long.reverseBytes(raFile.readLong()));
    }
    
    public double[] getValues(int off1, int off2) throws IOException { 
        if(off1 < start || off2 > end || off1 > off2) { 
            throw new IllegalArgumentException();            
        }
        
		System.out.println(">> " + off1 + "-" + off2 + " (" + start + "," + end +")");
        raFile.seek(sHeaderSize + ((long)(off1 - start) * sElementSize));
        int count = off2 - off1 + 1;
        double[] vals = new double[count];
		int i = 0;
		try { 
			for(i = 0; i < vals.length; i++) { 
				vals[i] = Double.longBitsToDouble(Long.reverseBytes(raFile.readLong()));
			}
		} catch(EOFException ef) { 
			System.out.println("** offset " + i + " produced an EOF.");	
			for(; i < vals.length; i++) { 
				vals[i] = 0.0;
			}
		}
        return vals;
    }
}
