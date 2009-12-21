/*
 * Author: tdanford
 * Date: Dec 19, 2008
 */
package edu.mit.csail.cgs.utils.models;

import java.io.*;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.json.*;

public interface ModelOutput<T extends Model> extends Closeable {

	public void writeModel(T m);
	public void flush();
	
	public static class LineWriter<T extends Model> implements ModelOutput<T> {
		
		private PrintStream ps;
		
		public LineWriter(OutputStream os) { 
			ps = new PrintStream(os);
		}

		public void close() {
			ps.close();
			ps = null;
		}

		public boolean isClosed() {
			return ps == null;
		}

		public void flush() {
			ps.flush();
		}

		public void writeModel(T m) {
			JSONObject obj = m.asJSON();
			ps.println(obj.toString());
		} 
		
	}
}
