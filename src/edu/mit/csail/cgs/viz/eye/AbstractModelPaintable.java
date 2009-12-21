/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.mit.csail.cgs.viz.eye;

import java.awt.*;
import java.util.*;

import javax.swing.Action;

import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.models.Model;
import edu.mit.csail.cgs.viz.paintable.AbstractPaintable;
import edu.mit.csail.cgs.viz.paintable.PaintableChangedListener;

public abstract class AbstractModelPaintable 
	extends AbstractPaintable implements ModelPaintable, ModelPaintablePropertyListener {
	
	private Map<String,ModelPaintableProperty> props;
	private Map<Point,Set<Model>> drawnPoints;
	private Map<Rectangle,Set<Model>> drawnRects;
	
	public AbstractModelPaintable() { 
		props = new HashMap<String,ModelPaintableProperty>();
		drawnPoints = null;
		drawnRects = null;
	}
	
	public void startDrawingPoints() { 
		drawnPoints = new HashMap<Point,Set<Model>>();
		drawnRects = new HashMap<Rectangle,Set<Model>>();
		System.out.println("Starting!!");
	}
	
	public void stopDrawingPoints() { 
		drawnPoints = null;
		drawnRects = null;
		System.out.println("Stopping!!");
	}
	
	protected void clearDrawnPoints() { 
		if(drawnPoints != null) { 
			drawnPoints.clear();
			drawnRects.clear();
		}
	}
	
	protected void drawPoint(Point p, Model m) { 
		if(drawnPoints != null) {
			if(!drawnPoints.containsKey(p)) { 
				drawnPoints.put(p, new HashSet<Model>());
			}
			drawnPoints.get(p).add(m);
		}
	}
	
	protected void drawRect(Rectangle r, Model m) { 
		if(drawnRects != null) { 
			if(!drawnRects.containsKey(r)){ 
				drawnRects.put(r, new HashSet<Model>());
			}
			drawnRects.get(r).add(m);
		}
	}
	
	public Collection<Pair<Rectangle,Model>> findDrawnRects(Point p) { 
		ArrayList<Pair<Rectangle,Model>> models = new ArrayList<Pair<Rectangle,Model>>();
		if(drawnRects != null) { 
			for(Rectangle r : drawnRects.keySet()) {
				if(r.contains(p)) { 
					for(Model m : drawnRects.get(r)) { 
						models.add(new Pair<Rectangle,Model>(r, m));
					}
				}
			}
		}
		return models;
	}
	
	public Collection<Model> findContainedDrawnPoints(Rectangle r) { 
		ArrayList<Model> models = new ArrayList<Model>();
		if(drawnPoints != null) { 
			for(Point p : drawnPoints.keySet()) { 
				if(r.contains(p)) { 
					models.addAll(drawnPoints.get(p));
				}
			}
		}
		return models;
	}
	
	public Pair<Point,Set<Model>> findNearestDrawnPoint(Point p) { 
		Set<Model> m = null;
		Double dist = null;
		Point pm = null;
		if(drawnPoints != null) { 
			for(Point pp : drawnPoints.keySet()) { 
				int dx = p.x-pp.x, dy = p.y-pp.y;
				double dpp = (dx*dx) + (dy*dy);
				if(dist == null || dpp < dist) { 
					dist = dpp;
					m = drawnPoints.get(pp);
					pm = pp;
				}
			}
		}
		return m != null ? new Pair<Point,Set<Model>>(pm, m) : null;
	}

	public abstract void paintItem(Graphics g, int x1, int y1, int x2, int y2);
	public abstract void addModel(Model m);
	public abstract void addModels(Iterator<? extends Model> itr);
	public abstract void clearModels();
	
	public void addValue(Object v) { 
		if(v instanceof Model) { 
			addModel((Model)v);
		} else { 
			throw new IllegalArgumentException(String.format("Unsupported value type: %s", 
					v.getClass().getSimpleName()));
		}
	}
	
	public <X> X getPropertyValue(String key) {
		if(!props.containsKey(key)) { 
			throw new IllegalArgumentException("Unknown value: " + key);
		}
		ModelPaintableProperty prop = getProperty(key);
		return (X)prop.getValue();
	}

	public <X> X getPropertyValue(String key, X def) { 
		if(props.containsKey(key)) { 
			return (X)getPropertyValue(key);
		} else { 
			return def;
		}
	}
	
	public ModelPaintable setProperty(String k, Object value) { 
		return setProperty(new PropertyValueWrapper(k, value));
	}

	public ModelPaintable setProperty(String key, ModelPaintableProperty p) {
		if(props.containsKey(key)) { 
			props.get(key).removeListener(this);
		}
		
		props.put(key, p);
		p.addListener(this);
		
		return this; 
	}
	
	public ModelPaintable setProperty(ModelPaintableProperty p) {
		return setProperty(p.getKey(), p);
	}
	
	protected void initProperty(ModelPaintableProperty p) {
		if(props.containsKey(p.getKey())) { 
			props.get(p.getKey()).removeListener(this);
		}
		
		props.put(p.getKey(), p);
		p.addListener(this);
	}
	
	public ModelPaintable synchronizeProperty(String k, ModelPaintable p) { 
		if(props.containsKey(k)) { 
			p.setProperty(props.get(k));
		}
		return this; 
	}
	
	public void synchronizeProperty(String k1, ModelPaintable p, String k2) { 
		if(props.containsKey(k1)) { 
			p.setProperty(k2, props.get(k1));
		}
	}
	
	public <T extends ModelPaintableProperty> T getProperty(String k) { 
		return (T)props.get(k);
	}
	
	public void propertyChanged(ModelPaintableProperty p) { 
		dispatchChangedEvent();
	}
}
