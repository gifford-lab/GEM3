/*
 * Author: tdanford
 * Date: Mar 21, 2009
 */
package edu.mit.csail.cgs.utils.models;

import java.lang.reflect.Array;

public class ArrayFieldAccessor<T extends Model> implements Accessor<T> {
	
	private Accessor<T> arrayAccessor;
	private int idx;
	
	public ArrayFieldAccessor(Accessor<T> acc, int i) { 
		arrayAccessor = acc;
		idx = i;
		if(!arrayAccessor.getType().isArray()) { 
			throw new IllegalArgumentException(String.format(
					"Type of %s isn't array.", arrayAccessor.getType().getSimpleName()));
		}
	}

	public ArrayFieldAccessor(Class<T> cls, String arrayName, int i) {
		this(new FieldAccessor<T>(cls, arrayName), i);
	}
	
	public int hashCode() { 
		int code = 17 + arrayAccessor.hashCode();
		code *= 37;
		code += idx; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof ArrayFieldAccessor)) { return false; }
		ArrayFieldAccessor f = (ArrayFieldAccessor)o;
		return f.arrayAccessor.equals(arrayAccessor) && idx == f.idx;
	}
	
	public String toString() { return getName(); }

	public Class getType() { 
		return arrayAccessor.getType().getComponentType();
	}
	
	public String getName() { 
		return String.format("%s[%d]", arrayAccessor.getName(), idx);
	}
	
	public String getBaseName() { 
		return arrayAccessor.getBaseName();
	}
	
	public void set(T object, Object value) { 
		Object array = arrayAccessor.get(object);
		Array.set(array, idx, value);
	}

	public Object get(T object) {
		Object array = arrayAccessor.get(object);
		return Array.get(array, idx);
	}
}
