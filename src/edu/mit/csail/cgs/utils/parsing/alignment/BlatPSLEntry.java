package edu.mit.csail.cgs.utils.parsing.alignment;

import java.io.*;
import java.util.*;

public class BlatPSLEntry extends AlignmentRecord {

	private int match, mismatch, repmatch;
	private int Ns, QgapCount, QgapBases, TgapCount, TgapBases;
	private char strand;
	private String Qname;
	private int Qsize, Qstart, Qend;
	private String Tname;
	private int Tsize, Tstart, Tend;
	private int blockCount;
	private Vector<Integer> blockSizes, qStarts, tStarts;
	
	public int getMatch() { return match; }
	public int getMismatch() { return mismatch; }
	public int getRepmatch() { return repmatch; }
	public int getNs() { return Ns; }
	public int getQgapCount() { return QgapCount; }
	public int getQgapBases() { return QgapBases; }
	public int getTgapCount() { return TgapCount; }
	public int getTgapBases() { return TgapBases; }
	public char getStrand() { return strand; }
	public String getQname() { return Qname; }
	public int getQsize() { return Qsize; }
	public int getQstart() { return Qstart; }
	public int getQend() { return Qend; }
	public String getTname() { return Tname; }
	public int getTsize() { return Tsize; }
	public int getTstart() { return Tstart; }
	public int getTend() { return Tend; }
	public int getBlockCount() { return blockCount; }
	public int getBlockSize(int i) { return blockSizes.get(i); }
	public int getQStart(int i) { return qStarts.get(i); }
	public int getTStart(int i) { return tStarts.get(i); }

	public int getTGapSize(int i) { 
		if(i < 0 || i >= blockSizes.size()-1) { 
			throw new IllegalArgumentException("gap: " + i + ", blocks: " + blockSizes.size()); 
		}
		return tStarts.get(i+1) - (tStarts.get(i) + blockSizes.get(i));
	}
	
	public int hashCode() { 
		int code = 17;
		code += match; code *= 37;
		code += mismatch; code *= 37;
		code += (strand == '+' ? 1 : 3); code *= 37;
		code += Qname.hashCode(); code *= 37;
		code += Qstart; code *= 37;
		code += Tname.hashCode(); code *= 37;
		code += Tstart; code *= 37;
		return code;
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof BlatPSLEntry)) { return false; }
		BlatPSLEntry e = (BlatPSLEntry)o;
		if(!Qname.equals(e.Qname)) { return false; }
		if(!Tname.equals(e.Tname)) { return false; }
		if(strand != e.strand) { return false; }
		if(Qstart != e.Qstart) { return false; }
		if(Qend != e.Qend) { return false; }
		if(Tstart != e.Tstart) { return false; }
		if(Tend != e.Tend) { return false; }
		if(match != e.match) { return false; }
		if(mismatch != e.mismatch) { return false; }
		if(repmatch != e.repmatch) { return false; }
		if(Ns != e.Ns) { return false; }
		if(QgapCount != e.QgapCount) { return false; }
		if(QgapBases != e.QgapBases) { return false; }
		if(TgapCount != e.TgapCount) { return false; }
		if(TgapBases != e.TgapBases) { return false; }
		if(blockCount != e.blockCount) { return false; }
		return true;
	}

	public BlatPSLEntry(String line) { 
		String[] array = line.split("\\t");
		
		match = Integer.parseInt(array[0]);
		mismatch = Integer.parseInt(array[1]);
		repmatch = Integer.parseInt(array[2]);
		
		Ns = Integer.parseInt(array[3]);
		QgapCount = Integer.parseInt(array[4]);
		QgapBases = Integer.parseInt(array[5]);
		TgapCount = Integer.parseInt(array[6]);
		TgapBases = Integer.parseInt(array[7]);
		
		strand = array[8].charAt(0);
		
		Qname = array[9];
		Qsize = Integer.parseInt(array[10]);
		Qstart = Integer.parseInt(array[11]);
		Qend = Integer.parseInt(array[12]);

		Tname = array[13];
		Tsize = Integer.parseInt(array[14]);
		Tstart = Integer.parseInt(array[15]);
		Tend = Integer.parseInt(array[16]);
		
		blockCount = Integer.parseInt(array[17]);
		
		blockSizes = new Vector<Integer>();
		qStarts = new Vector<Integer>();
		tStarts = new Vector<Integer>();
		
		String[] aa;
		aa = array[18].split(",");
		for(int i = 0; i < aa.length; i++) { 
			String aai = aa[i].trim();
			if(aai.length() > 0) {
				blockSizes.add(Integer.parseInt(aai));
			}
		}

		aa = array[19].split(",");
		for(int i = 0; i < aa.length; i++) { 
			String aai = aa[i].trim();
			if(aai.length() > 0) {
				qStarts.add(Integer.parseInt(aai));
			}
		}

		aa = array[20].split(",");
		for(int i = 0; i < aa.length; i++) { 
			String aai = aa[i].trim();
			if(aai.length() > 0) {
				tStarts.add(Integer.parseInt(aai));
			}
		}
}
}
