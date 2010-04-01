/*
 * Created on Sep 6, 2005
 */
package edu.mit.csail.cgs.utils.io.parsing.BIND;

import java.io.*;
import java.util.*;

/**
 * @author tdanford
 */
public class BindInteraction {
    
    private BindInteractor a, b;
    private int rgid, bid;
    private String AB;
    private Set<BindReference> refs;
    
    public BindInteraction(int rgid, int bid, BindInteractor _a, BindInteractor _b, String ab) { 
        a = _a;
        b = _b;
        this.rgid = rgid; this.bid = bid;
        AB = ab;
        refs = new HashSet<BindReference>();
    }
    
    public BindInteractor getA() { return a; }
    public BindInteractor getB() { return b; }
    public String getAB() { return AB; }
    public int getRGID() { return rgid; }
    public int getBID() { return bid; }
    public String getKey() { return rgid + "/" + bid; }
    public Collection<BindReference> getRefs() { return refs; }
    
    public void addReference(BindReference bref) { refs.add(bref); }

    public boolean involvesInteractor(BindInteractor bi) { 
        return a.equals(bi) || b.equals(bi);
    }
    
    public int hashCode() { 
        int code = 17;
        code += rgid; code *= 37;
        code += bid; code *= 37;
        return code;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof BindInteraction)) { return false; }
        BindInteraction bi = (BindInteraction)o;
        return rgid == bi.rgid && bid == bi.bid;
    }
    
    public String toString() { 
		StringBuilder sb = new StringBuilder();
		sb.append("[" + bid + "(" + rgid + ") ");
		sb.append(a.toString());
		sb.append(" --> ");
		sb.append(b.toString());
		sb.append("]");
		return sb.toString();
    }    
}
