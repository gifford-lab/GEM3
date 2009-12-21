/*
 * Author: tdanford
 * Date: Aug 13, 2008
 */
package edu.mit.csail.cgs.utils;

import java.util.*;

public abstract class Sorter {

	public static <Key extends Comparable,Value> void jointSort(Key[] keys, Value[] values) {
		if(keys.length != values.length) { 
			String msg = String.format("keys array length %d != values array length %d", 
					keys.length, values.length);
			throw new IllegalArgumentException(msg); 
		}
		Sortable[] s = new Sorter.Sortable[keys.length];
		for(int i = 0; i < s.length; i++) { s[i] = new Sortable<Key,Value>(keys[i], values[i]); }
		Arrays.sort(s);
		for(int i = 0; i < s.length; i++) {
			Sortable<Key,Value> si = (Sorter.Sortable<Key,Value>)s[i]; 
			keys[i] = si.key;
			values[i] = si.value;
		}
	}
	
	private static class Sortable<K extends Comparable,V> implements Comparable<Sortable<K,V>> {
		
		public K key;
		public V value;
		
		public Sortable(K k, V v) { key = k; value = v; }
		
		public int hashCode() { 
			int code = 17;
			code += key.hashCode(); code *= 37;
			return code;
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof Sortable)) { 
				return false; 
			}
			Sortable s = (Sortable)o;
			return s.key.equals(key) && s.value.equals(value);
		}

		public int compareTo(Sortable s) {
			return key.compareTo(s.key);
		} 
		
	}
}
