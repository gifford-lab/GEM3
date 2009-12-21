/*
 * Author: tdanford
 * Date: Nov 15, 2008
 */
package edu.mit.csail.cgs.utils.datastructures;

import java.util.*;

/**
 * An implementation of a Heap/priority-queue, which is basically derived from the 
 * pseudocode in: 
 * "Heapsort", Chapter 6 (pg. 127) of Cormen, Leiserson, Rivest, and Stein (ed. 2)
 * 
 * @author tdanford
 *
 * @param <Key>
 */
public class Heap<Key extends Comparable> {
	
	public static void main(String[] args) { 
		Integer[] array = new Integer[] { 5, 15, 23, 1, 2, 6, 21, 99 };
		Heap<Integer> heap = new Heap<Integer>(array, -1);
		
		while(heap.size() > 0) { 
			Integer s = heap.removeFirst();
			System.out.print(String.format("%d ", s));
		}
		System.out.println();
	}

	private int polarity; // == -1 or 1.  polarity == -1 -> min-heap, == 1 -> max-heap
	private Vector<KeyWrapper> array;
	
	public Heap(Heap<Key> h) { 
		array = new Vector<KeyWrapper>(h.array);
		polarity = h.polarity;
	}
	
	public Heap(Key[] unsorted, int dir) { 
		array = new Vector<KeyWrapper>();
		for(int i = 0; unsorted != null && i < unsorted.length; i++) { 
			array.add(new KeyWrapper(unsorted[i]));
		}
		
		polarity = dir;
		if(dir != -1 && dir != 1) { throw new IllegalArgumentException(String.valueOf(dir)); }
		for(int i = array.size()/2; i >= 0; i--) { 
			maxHeapify(i);
		}
	}
	
	public Heap() { 
		this(null, 1);
	}
	
	public Heap(int dir) { 
		this(null, dir);
	}
	
	public Heap(Key[] unsorted) { 
		this(unsorted, 1);
	}
	
	/*
	 * Public Methods
	 */
	
	public List<Key> asList() { 
		Heap<Key> h = new Heap<Key>(this);
		ArrayList<Key> keys = new ArrayList<Key>();
		
		while(h.size()>0) { 
			keys.add(h.removeFirst());
		}
		
		return keys;
	}
	
	public Key getFirst() { 
		return array.get(0).key;
	}
	
	public Key removeFirst() {
		if(array.isEmpty()) { 
			throw new IllegalArgumentException();
		}
		
		Key first = array.get(0).key;
		int last = array.size()-1;
		KeyWrapper w = array.remove(last);
		if(!array.isEmpty()) { 
			array.set(0, w);
			maxHeapify(0);
		}
		return first;
	}
	
	public int size() { return array.size(); }
	
	public void insert(Key k) { 
		KeyWrapper w = new KeyWrapper(-polarity);
		array.add(w);
		increaseKey(array.size()-1, k);
	}
	
	public void increase(Key k1, Key k2) { 
		for(int i = 0; i < array.size(); i++) {
			if(array.get(i).key.equals(k1)) { 
				increaseKey(i, k2);
				return;
			}
		}
		throw new IllegalArgumentException(k1.toString());
	}
	
	/*
	 * Helper Methods 
	 */
	
	private void increaseKey(int i, Key k) { 
		KeyWrapper old = array.get(i);
		KeyWrapper newWrapper = new KeyWrapper(k);
		if(newWrapper.compareTo(old) == -polarity) { 
			throw new IllegalArgumentException(
					String.format("Can't increase/decrease key %s over %s", 
							newWrapper.key.toString(), old.key.toString()));
		}
		
		array.set(i, newWrapper);
		while(i > 0 && array.get(parent(i)).compareTo(array.get(i)) == -polarity) {
			exchange(i, parent(i));
			i = parent(i);
		}
	}

	private int left(int i) { 
		return i << 1;
	}
	
	private int right(int i) { 
		return (i << 1) + 1;
	}
	
	private int parent(int i) { 
		return i >> 1;
	}
	
	private void maxHeapify(int i) { 
		int l = left(i), r = right(i);
		int extreme = i;
		if(l < array.size() && array.get(l).compareTo(array.get(i)) == polarity) { 
			extreme = l;
		}
		if(r < array.size() && array.get(r).compareTo(array.get(extreme)) == polarity) { 
			extreme = r;
		}
		
		if(extreme != i) {
			exchange(extreme, i);
			maxHeapify(extreme);
		}
	}
	
	private void exchange(int i1, int i2) { 
		KeyWrapper k1 = array.get(i1);
		array.set(i1, array.get(i2));
		array.set(i2, k1);
	}
	
	public static final int MIN = -1;
	public static final int MAX = 1;
	
	private class KeyWrapper implements Comparable<KeyWrapper> {
		
		public Key key;
		public int special;
		
		public KeyWrapper(Key k) { 
			key = k; 
			special = 0;
		}
		
		public KeyWrapper(int s) { 
			key = null; 
			special = s;
			if(special != MIN && special != MAX) { 
				throw new IllegalArgumentException(String.valueOf(special));
			}
		}
		
		public boolean equals(Object o) { 
			if(!(o instanceof Heap.KeyWrapper)) { return false; }
			KeyWrapper w = (KeyWrapper)o;
			return special == w.special && 
				(key == null || w.key == null ? key == w.key : key.equals(w.key));
		}
		
		public int hashCode() { 
			if(key != null) { 
				return key.hashCode(); 
			} else { 
				return (17 + special) * 37;
			}
		}
		
		public int compareTo(KeyWrapper kw) { 
			if(special == MIN) {
				if(kw.special==MIN) { 
					return 0; 
				} else { 
					return -1;
				}
			} else if (special == MAX) {
				if(kw.special==MAX) { 
					return 0; 
				} else { 
					return 1; 
				}
			} else { 
				if(kw.key != null) { 
					return key.compareTo(kw.key);
				} else { 
					return -kw.compareTo(this);
				}
			}
		}
	}
}