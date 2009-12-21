/*
 * Author: tdanford
 * Date: Mar 21, 2009
 */
package edu.mit.csail.cgs.utils.models;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;

public class FieldAccessor<T extends Model> implements Accessor<T> {
	
	private Field field;
	
	public FieldAccessor(Field f) { 
		field = f;
	}
	
	public FieldAccessor(Class cls, String fieldName) { 
		try {
			field = cls.getField(fieldName);
		} catch (NoSuchFieldException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(String.format("No field %s in class %s", fieldName, cls.getSimpleName()), e);
		}
		if((field.getModifiers() | Modifier.PUBLIC) == 0) { 
			throw new IllegalArgumentException(String.format("Class %s has field %s, " +
					"but the field isn't public.", cls.getSimpleName(), fieldName));
		}
	}
	
	public Class getType() { return field.getType(); }
	public String getName() { return field.getName(); }
	
	public String getBaseName() { 
		return field.getName();
	}
	
	public int hashCode() { 
		return field.hashCode();
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof FieldAccessor)) { return false; }
		FieldAccessor a = (FieldAccessor)o;
		return a.field.equals(field);
	}
	
	public String toString() { return getName(); }
	
	public void set(T object, Object value) { 
		try {
			field.set(object, value);
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(object.toString());
		}
	}

	public Object get(T object) {
		try {
			return field.get(object);
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(object.toString());
		}
	}

}
