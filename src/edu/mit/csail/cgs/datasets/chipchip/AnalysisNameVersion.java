package edu.mit.csail.cgs.datasets.chipchip;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class AnalysisNameVersion extends NameVersion implements Comparable<NameVersion> {
    
    public AnalysisNameVersion(DataInputStream dis) throws IOException { 
        super(dis);
    }
    
    public AnalysisNameVersion(AnalysisNameVersion env) { 
        super(env);
    }
    
    public AnalysisNameVersion(String n, String v) {
        super(n,v);
    }
    public int compareTo(AnalysisNameVersion nv) { 
        if(label != null && nv.label != null) { 
            return label.compareTo(nv.label);
        } 
        if(!name.equals(nv.name)) { return name.compareTo(nv.name); }
        if(!version.equals(nv.version)) { return version.compareTo(nv.version); }
        return 0;
    }
    
    public boolean equals(Object o) {
        if (! (o instanceof AnalysisNameVersion)) {
            return false;
        }
        AnalysisNameVersion other = (AnalysisNameVersion)o;
        return (other.name.equals(name) &&
                other.version.equals(version) &&
                ((label == null && other.label == null) ||
                 other.label.equals(label)));
    }
}
