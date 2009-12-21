package edu.mit.csail.cgs.viz.graphs;

import java.util.*;
import java.awt.*;

public class ObjectView implements View {
	
	protected Map<String,Object> options;

	public ObjectView() {
		options = new HashMap<String,Object>();
	}
	
	public ObjectView(ObjectView ov) { 
		options = new HashMap<String,Object>(ov.options);
	}
	
	public void paintView(Graphics2D g2) { 
		// do nothing.
	}

	public boolean containsOption(String k) { return options.containsKey(k); }
	public Object getOption(String k) { return options.get(k); }
	public void setOption(String k, Object v) { options.put(k, v); }
    public void clearOption(String k) { options.remove(k); }
	
	public String getName() { return (String)options.get("name"); }
	public void setName(String n) { options.put("name", n); }
	public String getID() { return (String)options.get("ID"); }
	public void setID(String n) { options.put("ID", n); }
	
	public int getX() { 
		if(containsOption("x")) { 
			return (Integer)options.get("x"); 
		} else { 
			return 0;
		}
	}
	
	public void setX(int x) { 
		options.put("x", x);
	}

	public int getY() { 
		if(containsOption("y")) { 
			return (Integer)options.get("y"); 
		} else { 
			return 0;
		}
	}
	
	public void setY(int y) { 
		options.put("y", y);
	}
	
	public Color getColor() { 
		if(containsOption("color")) { 
			return (Color)options.get("color");
		} else { 
			return Color.black;
		}
	}
	
	public void setColor(Color c) { 
		options.put("color", c);
	}
}
