package edu.mit.csail.cgs.viz.preferences;

import java.util.*;

public interface PreferencesModel {
	
	public Set<String> getKeys();
	public Object getValue(String key);
	public void setValue(String k, Object v);
	public void addListener(PreferencesListener pl);
	public void removeListener(PreferencesListener pl);
	
	public static class Default implements PreferencesModel {
		
		private LinkedHashMap<String,Object> values;
		private LinkedList<PreferencesListener> listeners;
		
		public Default() { 
			values = new LinkedHashMap<String,Object>();
			listeners = new LinkedList<PreferencesListener>();
		}
		
		public void fireUpdatedEvent(Object src) { 
			PreferencesEvent pe = new PreferencesEvent(src, this);
			for(PreferencesListener l : listeners) { 
				l.preferencesUpdated(pe);
			}
		}
		
		public void fireCanceledEvent(Object src) { 
			PreferencesEvent pe = new PreferencesEvent(src, this);
			for(PreferencesListener l : listeners) { 
				l.preferencesUpdateCanceled(pe);
			}			
		}
		
		public void setValue(String key, Object v) { 
			int sep = -1;
			if((sep = key.indexOf(":")) != -1) { 
				String start = key.substring(0, sep);
				String rest = key.substring(sep+1, key.length());
				
				if(!values.containsKey(start)) { 
					values.put(start, new PreferencesModel.Default());
				}
				
				PreferencesModel inner = (PreferencesModel)values.get(start);
				inner.setValue(rest, v);
			} else { 
				values.put(key, v);
			}
		}

		public void addListener(PreferencesListener pl) {
			listeners.addLast(pl);
		}

		public Set<String> getKeys() {
			LinkedHashSet<String> keys = new LinkedHashSet<String>();
			for(String k : values.keySet()) { 
				Object val = values.get(k);
				if(val instanceof PreferencesModel) {
					PreferencesModel inner = (PreferencesModel)val;
					for(String k2 : inner.getKeys()) { 
						keys.add(k + ":" + k2);
					}
				} else { 
					keys.add(k);
				}
			}
			return keys;
		}

		public Object getValue(String key) {
			int sep = -1;
			if((sep = key.indexOf(":")) != -1) { 
				String start = key.substring(0, sep);
				String rest = key.substring(sep+1, key.length());
				PreferencesModel inner = (PreferencesModel)values.get(start);
				return inner.getValue(rest);
			} else { 
				return values.get(key);
			}
		}

		public void removeListener(PreferencesListener pl) {
			listeners.remove(pl);
		} 
	}
}
