package edu.mit.csail.cgs.conservation;

import java.io.*;
import edu.mit.csail.cgs.utils.Saveable;

/**
 * @author tdanford
 * 
 * Describes an experiment which is part of the analysis.
 *
 */
public class ExptDescriptor implements Saveable, Comparable<ExptDescriptor> {
	
	private String name;

	public ExptDescriptor(String n) { 
		name = n; 
	}
	
	public ExptDescriptor(DataInputStream dis) throws IOException { 
		name = dis.readUTF();
	}
	
	public void save(DataOutputStream dos) throws IOException { 
		dos.writeUTF(name);
	}
	
	public String getName() { return name; }
	
	public String toString() { 
		return name;
	}
	
	public int hashCode() { 
		int code = 17;
		code += name.hashCode(); code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ExptDescriptor)) { return false; }
		ExptDescriptor ed = (ExptDescriptor)o;
		if(!name.equals(ed.name)) { return false; }
		return true;
	}

    /* (non-Javadoc)
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    public int compareTo(ExptDescriptor ed) {
        return name.compareTo(ed.name);
    }
}
