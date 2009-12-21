package edu.mit.csail.cgs.ewok.verbs.io;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;

import edu.mit.csail.cgs.ewok.verbs.Sink;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.CGSException;
import edu.mit.csail.cgs.datasets.alignments.parsing.BLAST8HitRegion;


public class BLAST8HitFileOutputStream implements Sink<BLAST8HitRegion>, Closeable {

	static Logger logger = Logger.getLogger(BLAST8HitFileOutputStream.class);

	private File outputFile = null;
	private FileOutputStream fos;
	private PrintStream ps;

	/**
	 * 
	 * @param filename
	 * @param organism
	 * @param genome
	 * @throws CGSException
	 */
	public BLAST8HitFileOutputStream(String filename) throws CGSException {
		this(new File(filename));
	}


	/**
	 * 
	 * @param f
	 * @param organism
	 * @param genome
	 * @throws CGSException
	 */
	public BLAST8HitFileOutputStream(File f) throws CGSException {		
		outputFile = f;
		try {
			fos = new FileOutputStream(outputFile);
			ps = new PrintStream(fos);
		}
		catch (IOException ioex) {
			throw new CGSException(ioex);
		}
	}
	
	
	/**
	 * 
	 */
	public void init() {}
	
	
	/**
	 * 
	 */
	public void consume(BLAST8HitRegion val) {
		ps.println(val.getBLAST8String());
	}

	
	/**
	 * 
	 * @param iter
	 */
	public void consume(Iterator<BLAST8HitRegion> iter) {
        while (iter.hasNext()) {
            consume(iter.next());
        }		
	}
	
	
	/**
	 * 
	 */
	public void finish() {
		close();
	}

	
	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.mit.csail.cgs.utils.Closeable#close()
	 */
	public void close() {
		if (isClosed()) {
			return;
		}
		try {
			ps.close();
			if (fos != null) {
				fos.close();
			}
			ps = null;
			fos = null;
		} catch (IOException ioex) {
			logger.error("Error closing buffered writer", ioex);

			// TODO throw a generic (e.g. CGS) exception

		}
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see edu.mit.csail.cgs.utils.Closeable#isClosed()
	 */
	public boolean isClosed() {
		return ps == null;
	}
}