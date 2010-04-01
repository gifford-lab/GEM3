/*
 * Author: tdanford
 * Date: Jun 20, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.awt.Color;
import java.util.*;
import java.util.regex.*;
import java.io.*;
import java.lang.reflect.*;

import edu.mit.csail.cgs.utils.json.JSONArray;
import edu.mit.csail.cgs.utils.json.JSONException;
import edu.mit.csail.cgs.utils.json.JSONObject;

public class Model {
	
	public Model() { 
	}
	
	public Model(JSONObject obj) { 
		setFromJSON(obj);
	}
	
	public Model(Map<String,Object> map) {
		setFromMap(map);
	}
	
	public Model(Model m) { 
		setFromModel(m);
	}

    public Collection<Field> getFields() {
        ArrayList<Field> output = new ArrayList<Field>();
		Field[] fs = getClass().getFields();
		for(int i = 0; i < fs.length; i++) {
			boolean isStatic = (fs[i].getModifiers() & Modifier.STATIC) != 0;
			if(!isStatic) {         
                output.add(fs[i]);
            }
        }
        return output;
    }
    
    public void save(File f) throws IOException { 
    	PrintStream ps = new PrintStream(new FileOutputStream(f));
    	ps.println(asJSON().toString());
    	ps.close();
    }
    
    public void load(File f) throws IOException { 
    	FileReader fr = new FileReader(f);
    	ModelInput rdr = new ModelInput.LineReader(getClass(), fr);
    	setFromModel(rdr.readModel());
    	fr.close();
    }
    
    public void updateModel() { 
    }

	public static boolean isSubclass(Class c1, Class c2) { 
		return c2.isAssignableFrom(c1);
	}
	
	public Model cloneModel() { 
		Class c = getClass(); 
		try {
			Model instance = (Model)c.newInstance();
			instance.setFromModel(this);
			return instance;
		} catch (InstantiationException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(c.getSimpleName());
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(c.getSimpleName());
		}
	}
	
	public void setFromMap(Map<String,Object> map) { 
		Class c = getClass();
		for(String k : map.keySet()) { 
			try {
				Field f = c.getField(k);
				f.set(this, map.get(k));
			} catch (NoSuchFieldException e) {
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		updateModel();
	}
	
	public void setFromModel(Model m) { 
		setFromModel(m, null);
	}
	
	public void setFromModel(Model m, String prefix) { 
		Field[] fs = getClass().getFields();
		Class mclass = m.getClass();
		
		Pattern prefixPattern =
			prefix != null ? 
			Pattern.compile(String.format("%s_(.*)", prefix)) : null;
		
		for(int i = 0; i < fs.length; i++) { 
			Class type = fs[i].getType();
			int modifier = fs[i].getModifiers();
			if((modifier & Modifier.STATIC) == 0 && 
			   (modifier & Modifier.FINAL) == 0) { 

				try {
					String fieldName = fs[i].getName();
					if(prefix != null) { 
						Matcher prefixMatch = prefixPattern.matcher(fieldName);
						if(prefixMatch.matches()) { 
							fieldName = prefixMatch.group(1);
						}
					}
					
					Field mfield = mclass.getField(fieldName);
					Class mtype = mfield.getType();
					if(type.isAssignableFrom(mtype)) { 
						fs[i].set(this, mfield.get(m));
					}

				} catch (NoSuchFieldException e) {
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
			}
		}
		
		updateModel();
	}
	
	public String toString() { return toString(-1); }
	
	public int hashCode() { 
		int code = 17;
		Field[] fields = getClass().getFields();
		for(int i = 0; i < fields.length; i++) {
			if((fields[i].getModifiers() & Modifier.STATIC) == 0 && 
			   (fields[i].getModifiers() & Modifier.PUBLIC) != 0) { 
				try {
					Object value = fields[i].get(this);
					Class type = fields[i].getType();
					if(value != null) { 
						if(Model.isSubclass(type, Integer.class)) { 
							code += (Integer)value;
						} else if(Model.isSubclass(type, Double.class)) { 
							long bits = Double.doubleToLongBits((Double)value);
							code += (int)(bits >> 32);
						} else if(!Model.isSubclass(type, Model.class)) {  
							code += value.hashCode();
						}
						code *= 37;
					}
				} catch (IllegalAccessException e) {
					//e.printStackTrace();
				}
			}
		}
		return code;
	}
	
	public boolean equals(Object o) { 
		Class mycls = getClass(), ocls = o.getClass();
		if(!mycls.equals(ocls)) { return false; }
		
		Field[] fields = getClass().getFields();
		for(int i = 0; i < fields.length; i++) {
			if((fields[i].getModifiers() & Modifier.STATIC) == 0 && 
			   (fields[i].getModifiers() & Modifier.PUBLIC) != 0) { 
				try {
					Object value = fields[i].get(this);
					Object ovalue = fields[i].get(o);
					
					Class type = fields[i].getType();
					if(Model.isSubclass(type, Model.class)) { 
						return value == ovalue;
					} else { 
						if(value != null && ovalue != null) { 
							if(!value.equals(ovalue)) { return false; }
						} else { 
							if(value != null || ovalue != null) { 
								return false; 
							}
						}
					}
				} catch (IllegalAccessException e) {
					//e.printStackTrace();
				}
			}
		}
		return true;
	}
	
	public static Object unjsonify(Class type, Object value) { 
		if(value == JSONObject.NULL) { 
			return null;
		} else if(value instanceof JSONObject) { 
			if(Model.isSubclass(type, Model.class)) { 
				Model m;
				try {
					m = (Model)(type.newInstance());
					m.setFromJSON((JSONObject)value);
					return m;
				} catch (InstantiationException e) {
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
				return null;
			} else { 
				return null;
			}
		} else if (value instanceof JSONArray) {
			if(type.isArray()) {
				
				Class arrayType = type.getComponentType();
				JSONArray jsonArray = (JSONArray)value;
				int length = jsonArray.length();
				Object array = Array.newInstance(arrayType, length);
				//System.err.print(String.format("\tArray(%d):", length));
				
				for(int i = 0; i < length; i++) { 
					try {
						Array.set(array, i, unjsonify(arrayType, jsonArray.get(i)));
						//System.err.print(" " + jsonArray.get(i));
					} catch (JSONException e) {
						e.printStackTrace();
						Array.set(array, i, null);
					}
				}
				
				//System.err.println();
				
				return array;
				
			} else { 
				System.err.println(type.toString() + " was not array-type.");
				return null;
			}
		} else if(Model.isSubclass(value.getClass(), type)) {
			return value;
        } else if (Model.isSubclass(type, java.lang.Number.class)) {
            if (isSubclass(type, java.lang.Double.class)) {
            	if(value instanceof String) { 
            		return Double.parseDouble((String)value);
            	} else { 
            		return new Double(((Number)value).doubleValue());
            	}
            } else if (isSubclass(type, java.lang.Integer.class)) {
            	if(value instanceof String) { 
            		return Integer.parseInt((String)value);
            	} else { 
            		return new Integer(((Number)value).intValue());
            	}
            } else if (isSubclass(type, java.lang.Float.class)) {
                return new Float(((Number)value).floatValue());
            } else if (isSubclass(type, java.lang.Long.class)) {
                return new Long(((Number)value).longValue());
            } else if (isSubclass(type, java.lang.Short.class)) {
                return new Short(((Number)value).shortValue());
            } else if (isSubclass(type, java.lang.Byte.class)) {
                return new Byte(((Number)value).byteValue());
            }           
		} else if (type.getName().equals("java.awt.Color")) {
            String[] pieces = ((String)value).split(",");
            return new Color(Integer.parseInt(pieces[0]),
                             Integer.parseInt(pieces[1]),
                             Integer.parseInt(pieces[2]),
                             Integer.parseInt(pieces[3]));
        }

		throw new RuntimeException("Couldn't map " + value + " to " + type);
        //		return null;
	}
	
	public static Object jsonify(Object value) throws JSONException {
		if(value == null) { 
			return JSONObject.NULL;
		} 
		
		Class vclass = value.getClass();

		if(vclass.isArray()) {
			JSONArray array = new JSONArray();
			int length = Array.getLength(value);
			
			for(int i = 0; i < length; i++) { 
				array.put(i, jsonify(Array.get(value, i)));
			}
			
			return array;
			
		} else if(Model.isSubclass(vclass, Boolean.class)) { 
			return value;
		} else if(Model.isSubclass(vclass, String.class)) { 
			return value;
		} else if (Model.isSubclass(vclass, Integer.class)) { 
			return value;
		} else if (Model.isSubclass(vclass, Double.class)) { 
			return value;
		} else if (Model.isSubclass(vclass, Long.class)) { 
			return value;
		} else if (Model.isSubclass(vclass, Model.class)) { 
			return ((Model)value).asJSON();
		} else if (Model.isSubclass(vclass, Color.class)) {
            Color c = (Color)value;
            return String.format("%d,%d,%d,%d", c.getRed(), c.getGreen(), c.getBlue(), c.getAlpha());
        }


		
		return null;
	}

	public JSONObject asJSON() { 
		JSONObject obj = new JSONObject();

		Field[] fs = getClass().getFields();
		for(int i = 0; i < fs.length; i++) {
			String fname = fs[i].getName();
			boolean isStatic = (fs[i].getModifiers() & Modifier.STATIC) != 0;
			if(!isStatic) { 
				try {
					Object fvalue = fs[i].get(this);
					obj.put(fname, jsonify(fvalue));

				} catch (IllegalAccessException e) {
					e.printStackTrace(System.err);
				} catch (edu.mit.csail.cgs.utils.json.JSONException e) {
					System.err.println(String.format("Model Class: %s", getClass().getName()));
					System.err.println(String.format("Field Name: %s (%s)", fs[i].getName(), fs[i].getType().getName()));
					System.err.println(String.format("Value: %s", toString()));
					e.printStackTrace(System.err);
				}
			}
		}

		return obj;
	}
	
	public void setFromJSON(JSONObject obj) { 
		//System.err.println(String.format("setFromJSON(%s)", asJSON()));
		
		Field[] fs = getClass().getFields();
		for(int i = 0; i < fs.length; i++) {
			boolean isStatic = (fs[i].getModifiers() & Modifier.STATIC) != 0;
			if(!isStatic) { 
				try {
					String fname = fs[i].getName();
					if(obj.has(fname)) { 
						Object value = obj.get(fname);
                        Class type = fs[i].getType();
                        Object unjsoned = unjsonify(fs[i].getType(), value);                        
                        fs[i].set(this, unjsoned);
						//System.err.println(String.format("\t%s <- %s (%s)", fs[i].getName(), unjsoned, fs[i].get(this)));
					} else { 
                        /* this used to set object fields to null if they weren't present
                           in the json representation.  I've changed it to ignore
                           those fields.  The old behavior was a problem if you were reading
                           property objects back from disk and the file was missing fields
                           that were recently added to an object.  Those fields would
                           be set to null even if the object's constructor had initialized
                           them to something else.
                        */
                        //						fs[i].set(this, null);
					}
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				} catch (JSONException e) {
					e.printStackTrace();
				}
			}
		}
		
		//System.err.println(String.format("\t->%s", asJSON()));
		updateModel();
	}
	
	public String toString(int maxLineLength) { 
		Field[] fs = getClass().getFields();
		StringBuilder sb = new StringBuilder();
		int linelength = 0;

		for(int i = 0; i < fs.length; i++) {
			try {
				String str = "";
				if(fs[i].getType().isArray()) { 
					Object array = fs[i].get(this);
					if(array != null) { 
						int len = Array.getLength(array);
						StringBuilder ab = new StringBuilder();
						ab.append("[");

						for(int j = 0; j < len; j++) { 
							Object value = Array.get(array, j);
							String valueStr = 
								value instanceof Double ? 
										String.format("%.3f", (Double)value) : String.valueOf(value);
										ab.append(String.format(" %s", valueStr));
						}
						ab.append(" ]");
						str = String.format("%s:%s", fs[i].getName(), ab.toString());
					} else { 
						str = String.format("%s:%s", fs[i].getName(), null);
					}

				} else {
					Object value = fs[i].get(this);
					String valueStr = 
						value instanceof Double ? 
								String.format("%.3f", (Double)value) : String.valueOf(value);
								str = String.format("%s:%s", fs[i].getName(), valueStr);
				}
				
				if(maxLineLength > 0 && linelength + str.length() > maxLineLength) { 
					sb.append("\n");
					linelength = 0;
				}

				sb.append(str);
				sb.append(" ");
				linelength += str.length() + 1;
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}

		return sb.toString();
	}

	public static class FieldComparator<T extends Model> implements Comparator<T> { 
		private Field[] fields;
		
		public FieldComparator(Class<T> cls, String... fnames) { 
			ModelFieldAnalysis<T> analysis = new ModelFieldAnalysis<T>(cls);
			Vector<Field> fs = new Vector<Field>();
			for(int i = 0; i < fnames.length; i++) { 
				Field f = analysis.findField(fnames[i]);
				if(f != null){ 
					fs.add(f);
				}
			}
			fields = fs.toArray(new Field[0]);
		}
		
		public FieldComparator(Field... fs) { 
			fields = fs.clone();
		}
		
		public int compare(T t1, T t2) {
			for(int i = 0; i < fields.length; i++) { 
				if(fields[i] != null) { 
					Class type = fields[i].getType();
					if(Model.isSubclass(type, Comparable.class)) { 
						try {
							Comparable v1 = (Comparable)fields[i].get(t1);
							Comparable v2 = (Comparable)fields[i].get(t2);
							int c = v1.compareTo(v2);
							if(c != 0) { return c; }
						} catch (IllegalAccessException e) {
							// do nothing!
						}
					}
				}
			}
			return 0;
		}
	}

	public static <T extends Model> T loadFromFile(Class<T> c, File f) throws IOException {
		FileReader reader = new FileReader(f);
		ModelInput<T> input = new ModelInput.LineReader<T>(c, reader);
		T m = input.readModel();
		reader.close();
		return m;
	}
	
	public static <T extends Model> void writeToFile(T model, File f) throws IOException { 
		PrintStream ps = new PrintStream(new FileOutputStream(f));
		ps.println(model.asJSON().toString());
		ps.close();
	}
}
