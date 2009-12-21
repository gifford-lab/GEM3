/*
 * Author: tdanford
 * Date: Nov 14, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.util.*;

import edu.mit.csail.cgs.utils.models.Model;

public class ContainerModelPaintable extends AbstractModelPaintable {
	
	protected ArrayList<ModelPaintable> innerPaintables;
	
	public ContainerModelPaintable() { 
		innerPaintables = new ArrayList<ModelPaintable>();
	}
	
	public ContainerModelPaintable(ModelPaintable... pts) { 
		this();
		
		for(int i = 0; i < pts.length; i++) { 
			addModelPaintable(pts[i]);
		}
	}
	
	public ContainerModelPaintable(Collection<ModelPaintable> pts) { 
		this();
		
		for(ModelPaintable p : pts) { 
			addModelPaintable(p);
		}
	}
	
	public void addModelPaintable(ModelPaintable p) {
		if(!innerPaintables.contains(p)) { 
			p.addPaintableChangedListener(this);
		}
		innerPaintables.add(p);
	}
	
	public void removeModelPaintable(ModelPaintable p) { 
		innerPaintables.remove(p);
		if(!innerPaintables.contains(p)) { 
			p.removePaintableChangedListener(this);
		}
	}

	public void addModel(Model m) {
		setEventPassthrough(false);
		for(ModelPaintable p : innerPaintables) { 
			p.addModel(m);
		}
		setEventPassthrough(true);
		dispatchChangedEvent();
	}

	public void addModels(Iterator<? extends Model> itr) {
		setEventPassthrough(false);
		while(itr.hasNext()) { 
			Model m = itr.next();
			for(ModelPaintable p : innerPaintables) { 
				p.addModel(m);
			}
		}
		setEventPassthrough(true);
		dispatchChangedEvent();
	}

	public void clearModels() {
		setEventPassthrough(false);
		for(ModelPaintable p : innerPaintables) { 
			p.clearModels();
		}
		setEventPassthrough(true);
		dispatchChangedEvent();		
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		// does nothing, by default.
	}
}
