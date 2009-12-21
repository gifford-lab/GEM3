/*
 * Author: tdanford
 * Date: Aug 20, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.util.*;
import java.awt.Color;
import java.lang.reflect.*;


public class ModelFieldAnalysis<T extends Model> {
	
	protected Class<? extends T> modelClass;
	protected Vector<Field> fields;
	
	protected Map<String,Field> staticFields;

	/**
	 * @param cls A Class object that represents a java class which is a subclass of Model. 
	 */
	public ModelFieldAnalysis(Class<? extends T> cls) { 
		modelClass = cls;
		fields = new Vector<Field>();
		staticFields = new HashMap<String,Field>();
		
		if(Model.isSubclass(modelClass, Model.class)) { 
			Field[] fa = modelClass.getFields();
			for(int i = 0; i < fa.length; i++) { 
				Class type = fa[i].getType();
				if((fa[i].getModifiers() & Modifier.STATIC) == 0) { 
					fields.add(fa[i]);
				} else { 
					staticFields.put(fa[i].getName(), fa[i]);
				}
			}
		} else {
			String msg = String.format("%s is not a subclass of Model", modelClass.getName());
			throw new IllegalArgumentException(msg);
		}
	}
	
	public Object get(String compositeFieldName, Object value) { 
		String[] fieldNameArray = compositeFieldName.split("\\.");
		
		for(int i = 0; i < fieldNameArray.length; i++) { 
			Class valueClass = value.getClass();
			String fieldName = fieldNameArray[i];
			if(!Model.isSubclass(valueClass, Model.class)) { return null; }
			try {
				Field field = valueClass.getField(fieldName);
				int modifiers = field.getModifiers();
				if((modifiers & Modifier.STATIC) != 0) { return null; }
				
				if(field == null) { return null; }
				if(value == null) { return null; }
				
				value = field.get(value);
				
			} catch (IllegalAccessException e) {
				return null;
			} catch (NoSuchFieldException e) {
				return null;
			}
		}
		
		return value;
	}
	
	public Class<? extends T> getModelClass() { return modelClass; }
	
	public Vector<Field> getFields() { return fields; }
	
	public boolean getStaticSwitch(String name, boolean defaultValue) { 
		if(!staticFields.containsKey(name)) { return defaultValue; }
		if(!Model.isSubclass(staticFields.get(name).getType(), Boolean.class)) { return defaultValue; }
		try {
			Boolean val = (Boolean) staticFields.get(name).get(modelClass);
			return val;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(String.format("Can't access field %s in class %s",
					name, modelClass.getSimpleName()));
		}
	}
	
	public Field[] findFields(String... fieldNames) { 
		Field[] array = new Field[fieldNames.length];
		for(int i = 0; i < fieldNames.length; i++) { 
			array[i] = findField(fieldNames[i]);
		}
		return array;
	}
	
	public Field findField(String fieldName) { 
		try {
			return modelClass.getField(fieldName);
		} catch (SecurityException e) {
			return null;
		} catch (NoSuchFieldException e) {
			return null;
		}
	}
	
	public Field findStaticField(String fieldName) { 
		return staticFields.get(fieldName); 
	}
	
	public Vector<String> getFieldNames() { 
		Vector<String> fns = new Vector<String>();
		for(Field f : fields) { 
			fns.add(f.getName());
		}
		return fns;
	}

	public Field findTypedField(String name, Class c) {
		for(int i = 0; i < fields.size(); i++) { 
			if((fields.get(i).getModifiers() & Modifier.STATIC) == 0) { 
				Class type = fields.get(i).getType();
				if(fields.get(i).getName().equals(name) && Model.isSubclass(type, c)) { 
					return fields.get(i);
				}
			}
		}
		return null;
	}

	public Vector<Field> findTypedFields(Class c) { 
		Vector<Field> fs = new Vector<Field>();
		for(int i = 0; i < fields.size(); i++) { 
			if((fields.get(i).getModifiers() & Modifier.STATIC) == 0) { 
				Class type = fields.get(i).getType();
				if(Model.isSubclass(type, c)) { 
					fs.add(fields.get(i));
				}
			}
		}
		return fs;		
	}

    public Vector<Field> findArrayFields() {
        Vector<Field> fs = new Vector<Field>();
        for(int i = 0; i < fields.size(); i++) { 
			if((fields.get(i).getModifiers() & Modifier.STATIC) == 0) { 
				Class type = fields.get(i).getType();
                if (type.isArray()) {
                    fs.add(fields.get(i));
                }
            }
        }
        return fs;
    }
    
    public Vector<Field> findTypedArrayFields(Class c) { 
        Vector<Field> fs = new Vector<Field>();
        for(int i = 0; i < fields.size(); i++) { 
			if((fields.get(i).getModifiers() & Modifier.STATIC) == 0) { 
				Class type = fields.get(i).getType();
                if (type.isArray() && Model.isSubclass(type.getComponentType(), c)) { 
                    fs.add(fields.get(i));
                }
            }
        }
        return fs;    	
    }
	
	public Vector<Field> foreignKeyFields() {
		return findTypedFields(Model.class);
	}
	
	public Vector<Field> nonForeignKeyFields() { 
		Vector<Field> fs = new Vector<Field>();
		for(int i = 0; i < fields.size(); i++) { 
			if((fields.get(i).getModifiers() & Modifier.STATIC) == 0) { 
				Class type = fields.get(i).getType();
				if(!Model.isSubclass(type, Model.class)) { 
					fs.add(fields.get(i));
				}
			}
		}
		return fs;
	}
}
