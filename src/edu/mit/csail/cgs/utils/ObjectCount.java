package edu.mit.csail.cgs.utils;

import java.lang.*;
import java.io.*;
import java.util.*;

public class ObjectCount implements Comparable { 

	private String fToken;
	private int fCount;

	public ObjectCount(String token, int c) { 
		fToken = token;
		fCount = c;
	}

	public ObjectCount(String token) { 
		fToken = token;
		fCount = 1;
	}

	public ObjectCount(ObjectCount oc) { 
		fToken = oc.fToken;
		fCount = oc.fCount;
	}

	public int getCount() { return fCount; }
	public void changeCount(int c) { fCount += c; }
	public String getToken() { return fToken; }

	// necessarily doesn't include the count, since this could 
	// change, which would require re-hashing the count in all
	// containing Counted objects, a big hassle.
	public int hashCode() { 
		int code = 17;
		code += fToken.hashCode(); code *= 37;
		return code;
	}

	public int compareTo(Object o) { 
		ObjectCount oc = (ObjectCount)o;
		if(fCount < oc.fCount) { return -1; }
		if(fCount > oc.fCount) { return 1; }
		return fToken.compareTo(oc.fToken);
	}

	public boolean equals(Object o) { 
		if(!(o instanceof ObjectCount)) { return false; }
		ObjectCount oc = (ObjectCount)o;
		if(!fToken.equals(oc.fToken)) { return false; }
		if(fCount != oc.fCount) { return false; }
		return true;
	}

	public String toString() { return "[" + fToken + "] : " + fCount; }
}

