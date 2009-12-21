package edu.mit.csail.cgs.datasets.chipchip;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class NameVersion implements Comparable<NameVersion> {
    public String name, version, label;
    
    public NameVersion(DataInputStream dis) throws IOException { 
        String n = dis.readUTF();
        String v = dis.readUTF();
        String l = dis.readUTF();
        if(l.length() == 0) { l = null; }        
        init(n,v,l);
    }
    
    public NameVersion(NameVersion nv) { 
        init(nv.name,nv.version,nv.label);
    }
    
    public NameVersion(String n, String v) {
        init(n,v,null);
    }
    public NameVersion(String n, String v, String l) {
        init(n,v,l);
    }
    public String getName() {return name;}
    public String getVersion() {return version;}
    public String getLabel() {return label;}
    public void init(String n, String v, String l) {
        if (n == null) {
            throw new NullPointerException("Name can't be null in a NameVersion");
        }
        if (v == null) {
            throw new NullPointerException("Version can't be null in a NameVersion");
        }
        name = n;
        version = v;
        label = l;
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        dos.writeUTF(name);
        dos.writeUTF(version);
        if(label==null) { 
            dos.writeUTF("");
        } else { 
            dos.writeUTF(label);
        }
    }

    public void setLabel(String l) {
        label = l;
    }
    
    public String toString() {
        return (label == null?"":label + ": ") + name + "(" + version + ")";
    }
    
    public int hashCode() { 
    	int code = 17;
    	code += name.hashCode(); code *= 37;
    	code += version.hashCode(); code *= 37;
    	if(label != null) { code += label.hashCode(); code *= 37; }
    	return code;
    }
    
    public int compareTo(NameVersion nv) { 
        if(!name.equals(nv.name)) { return name.compareTo(nv.name); }
        if(!version.equals(nv.version)) { return version.compareTo(nv.version); }
        if(label != null && nv.label != null) { 
            return label.compareTo(nv.label);
        } 
        return 0;
    }
    
    public boolean equals(Object o) {
        if (! (o instanceof NameVersion)) {
            return false;
        }
        NameVersion other = (NameVersion)o;
        if(!other.name.equals(name)) { return false; }
        if(!other.version.equals(version)) { return false; }
        if(label != null && other.label != null) {
            return label.equals(other.label);
        }
        return (label == null) == (other.label == null);
    }
}
