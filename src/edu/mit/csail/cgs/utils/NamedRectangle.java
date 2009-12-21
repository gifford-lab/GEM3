/*
 * Author: tdanford
 * Date: Sep 29, 2008
 */
package edu.mit.csail.cgs.utils;

import java.awt.Rectangle;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.utils.models.Modeleable;

public class NamedRectangle extends Rectangle implements Modeleable {
	
	private String name;
	
	public NamedRectangle(NamedRectangleModel m) { 
		super(m.x, m.y, m.width, m.height);
		name = m.name;
	}
	
	public NamedRectangle(int x, int y, int w, int h, String n) { 
		super(x, y, w, h);
		name = n;
	}
	
	public NamedRectangle(Rectangle r, String n) { 
		super(r);
		name = n;
	}
	
	public NamedRectangle(DataInputStream dis) throws IOException { 
		super(dis.readInt(), dis.readInt(), dis.readInt(), dis.readInt());
		name = dis.readUTF();
	}
	
	public String getName() { return name; }
	
	public int area() { 
		return width*height;
	}
	
	public int hashCode() { return name.hashCode(); }
	public boolean equals(Object o) { 
		if(!(o instanceof NamedRectangle)) { return false; }
		NamedRectangle r = (NamedRectangle)o;
		if(!name.equals(r.name)) { return false; }
		if(!super.equals(r)) { return false; }
		return true;
	}
	
	public void save(DataOutputStream dos) throws IOException { 
		dos.writeInt(x);
		dos.writeInt(y);
		dos.writeInt(width);
		dos.writeInt(height);
		dos.writeUTF(name);
	}
	
	public Class getModelClass() { return NamedRectangleModel.class; }
	
	public Model asModel() { 
		return new NamedRectangleModel(x, y, width, height, name);
	}
	
	public static class NamedRectangleModel extends Model { 
		public Integer x, y, width, height;
		public String name;
		
		public NamedRectangleModel(Integer x, Integer y, Integer w, Integer h, String n) { 
			this.x = x;
			this.y = y;
			width = w;
			height = h;
			name = n;
		}
	}
}