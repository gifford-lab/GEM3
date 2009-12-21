/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package edu.mit.csail.cgs.utils.models;

import java.io.*;
import java.util.*;

public class FileTimer implements Timer, edu.mit.csail.cgs.utils.Closeable {
	
	private PrintStream ps;
	
	public FileTimer(File f) throws IOException { 
		ps = new PrintStream(new FileOutputStream(f));
	}

	public FileTimer(File f, boolean append) throws IOException { 
		ps = new PrintStream(new FileOutputStream(f, append));
	}

	public void addTiming(Timing t) {
		ps.println(t.asJSON().toString());
	}

	public void addTimings(Iterator<Timing> ts) {
		while(ts.hasNext()) { 
			addTiming(ts.next());
		}
	}

	public void close() {
		ps.close();
		ps = null;
	}

	public boolean isClosed() {
		return ps == null;
	}
}
