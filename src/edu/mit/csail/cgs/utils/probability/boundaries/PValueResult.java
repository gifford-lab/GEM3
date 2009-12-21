/*
 * Created on Feb 14, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;
import java.io.*;
import java.text.*;

import edu.mit.csail.cgs.utils.*;

/**
 * @author tdanford
 */
public class PValueResult implements Saveable {
    
    private static NumberFormat nf;
    
    static { 
        nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(3);
    }
    
    public static PValueResult getResult(PValuator pv, BoundaryDataset ds) { 
        if(!ds.hasBoundary()) { ds.findBoundary(); }
        PValueResult pvr = pv.logPValue(ds);
        return pvr;
    }
    
    public static Collection<PValueResult> getResults(PValuator pv,
            Collection<BoundaryDataset> dsList) {
        
        LinkedList<PValueResult> pvList = new LinkedList<PValueResult>();
        for(BoundaryDataset ds : dsList) { 
            pvList.addLast(getResult(pv, ds));
        }
        return pvList;
    }
    
    private String descriptor;
    private int direction, error;
    private double logPValue;

    public PValueResult(String desc, int e, int dir, double lpv) {
        descriptor = desc;
        error = e;
        direction = dir;
        logPValue = lpv;
    }
    
    public PValueResult(DataInputStream dis) throws IOException { 
        descriptor = dis.readUTF();
        direction = dis.readInt();
        error = dis.readInt();
        logPValue = dis.readDouble();
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        dos.writeUTF(descriptor);
        dos.writeInt(direction);
        dos.writeInt(error);
        dos.writeDouble(logPValue);
    }
    
    public String getDescriptor() { return descriptor; }
    public int getError() { return error; }
    public int getDirection() { return direction; }
    public double getLogPValue() { return logPValue; }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        if(direction == 1) { 
            sb.append("+");
        } else { 
            sb.append("-");
        }
        sb.append(" ");
        sb.append(String.valueOf(error)); 
        sb.append("e ");
        sb.append(nf.format(logPValue) + " (" + nf.format(Math.exp(logPValue)) + ")");
        sb.append(" " + descriptor);
        return sb.toString();
    }
}
