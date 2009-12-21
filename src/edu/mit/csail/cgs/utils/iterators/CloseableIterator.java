/*
 * Author: tdanford
 * Date: Aug 8, 2008
 */
package edu.mit.csail.cgs.utils.iterators;

import java.util.Iterator;

import edu.mit.csail.cgs.utils.Closeable;

public interface CloseableIterator<X> extends Iterator<X>, Closeable {
	
	public static class Wrapper<Y> implements CloseableIterator<Y> {
		
		private Iterator<Y> itr;
		
		public Wrapper(Iterator<Y> i) { itr = i; }

		public boolean hasNext() {
			return itr.hasNext();
		}

		public Y next() {
			return itr.next();
		}

		public void remove() {
			itr.remove();
		}

		public void close() {
			if(itr instanceof Closeable) { 
				((Closeable)itr).close();
			}
			itr = null;
		}

		public boolean isClosed() {
			return itr == null;
		} 
		
	}
}
