package edu.mit.csail.cgs.datasets.binding;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.datasets.general.Typed;

import java.io.*;

public class BindingEvent extends Region implements Typed {

    private double size, conf;
    private String type;
    
    public BindingEvent(BindingEvent evt) { 
        super(evt.getGenome(), evt.getChrom(), evt.getStart(), evt.getEnd());
        size = evt.size;
        conf = evt.conf;
        type = evt.type;
    }
    
    public BindingEvent(BindingExtent ext) { 
    	super(ext.getGenome(), ext.getChrom(), ext.getExtentStart(), ext.getExtentEnd());
    	size = ext.getSize();
    	conf = ext.getConf();
    	type = ext.getType();
    }

    public BindingEvent(Genome g, String c, int start, int end, double size, double conf, String type) {
        super(g,c,start,end);
        this.size = size;
        this.conf = conf;
        this.type = type;
    }
    
    public BindingEvent(Region r, double s, double c, String t) { 
    	super(r);
    	this.size = s;
    	this.conf = c;
    	this.type = t;
    }
    
    public BindingEvent(Genome g, DataInputStream dis) throws IOException { 
        super(g, dis);
        size = dis.readDouble();
        conf = dis.readDouble();
        type = dis.readUTF();
    }
    
    public void save(DataOutputStream dos) throws IOException { 
        super.save(dos);
        dos.writeDouble(size);
        dos.writeDouble(conf);
        dos.writeUTF(type);
    }

    public double getSize() { return size; }
    public double getConf() { return conf; }
    public String getType() { return type; }

    public BindingEvent combine(BindingEvent e) { 
        if(!type.equals(e.type)) { 
            throw new IllegalArgumentException(type + " != " + e.type); 
        }
        
        if(!getChrom().equals(e.getChrom())) { 
            throw new IllegalArgumentException(getChrom() + " != " + e.getChrom()); 
        }
        
        int ns = Math.min(getStart(), e.getStart());
        int ne = Math.max(getEnd(), e.getEnd());
        double nsize = Math.max(size, e.size);
        double nconf = Math.max(conf, e.conf);
        
        BindingEvent newEvent = new BindingEvent(getGenome(), getChrom(), ns, ne, nsize, nconf, type);
        return newEvent;
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder(super.toString());
        sb.append(" [" + type + "]");
        sb.append(" size: " + size);
        sb.append(", conf: " + conf);
        return sb.toString();
    }
    
    public boolean equals(Object o) { 
        if(!super.equals(o)) { return false; }
        if(!(o instanceof BindingEvent)) { return false; }
        BindingEvent e = (BindingEvent)o;
        if(!type.equals(e.type)) { return false; }
        if(size != e.size) { return false; }
        if(conf != e.conf) { return false; }
        return true;
    }
    
    public int hashCode() { 
        int code = super.hashCode();
        
        code += type.hashCode(); code *= 37;
        long bits = Double.doubleToLongBits(size);
        code += (int)(bits >> 32); code *= 37;
        bits = Double.doubleToLongBits(conf);
        code += (int)(bits >> 32); code *= 37;
        
        return code;
    }
}
