package edu.mit.csail.cgs.datasets.chipchip;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class ExptNameVersion extends NameVersion implements Comparable<NameVersion> {
    public String replicate;
    
    public ExptNameVersion(DataInputStream dis) throws IOException { 
        super(dis);
        String r = dis.readUTF();
        if (r.length() == 0) {replicate = null;} else {replicate = r;}
    }
    
    public ExptNameVersion(ExptNameVersion env) { 
        super(env);
        replicate = env.replicate;
    }
    
    public ExptNameVersion(String n, String v) {
        super(n,v);
        replicate = null;
    }
    public ExptNameVersion(String n, String v, String r) {
        super(n,v);
        replicate = r;
    }
    public ExptNameVersion(String pieces[]) {
        super(pieces[0],pieces[1]);
        if (pieces.length == 3) {
            replicate = pieces[2];
        } 
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        super.save(dos);
        if(replicate==null) { 
            dos.writeUTF("");
        } else { 
            dos.writeUTF(replicate);
        }
    }
    
    public String getReplicate() {
        return replicate;
    }
    protected void setReplicate(String r) {replicate =r ;}
    public void clearReplicate() {replicate = null;}
   
    public String toString() {
        if (replicate != null) {
            return (label == null?"":label + ": ") + name + "(" + version + "," + 
                replicate + ")";
        } else {
            return super.toString();
        }
    }
    
    public int hashCode() { 
        int code = super.hashCode() * 37;
        if (replicate != null) {code += replicate.hashCode();}
        return code;
    }
    
    public int compareTo(NameVersion nv) { 
        int c = super.compareTo(nv);
        if (nv instanceof ExptNameVersion) {
            ExptNameVersion env = (ExptNameVersion) nv;
            if (c == 0) {
                if (replicate != null && env.replicate != null) {
                    c = replicate.compareTo(env.replicate);
                } else {
                    c = replicate == null ? -1 : 1;
                }
            } 
        }
        return c;
    }
    
    public boolean equals(Object o) {
        if (! (o instanceof ExptNameVersion)) {
            return false;
        }
        ExptNameVersion other = (ExptNameVersion)o;
        if(!super.equals(o)) { return false; }
        if(replicate != null && other.replicate != null) { 
            return replicate.equals(other.replicate);
        }
        return (replicate == null) == (other.replicate == null);
    }
}
