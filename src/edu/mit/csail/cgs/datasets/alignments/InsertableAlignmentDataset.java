package edu.mit.csail.cgs.datasets.alignments;

import java.util.*;
import java.io.*;

import edu.mit.csail.cgs.datasets.general.*;
import edu.mit.csail.cgs.datasets.species.*;

public class InsertableAlignmentDataset {

	public InsertableAlignmentVersion version;
	public Vector<InsertableAlignment> alignments;
	public Vector<InsertableAlignBlock> blocks;
	
	public InsertableAlignmentDataset(InsertableAlignmentVersion v) { 
		version = v;
		alignments = new Vector<InsertableAlignment>();
		blocks = new Vector<InsertableAlignBlock>();
	}
	
	public boolean check() { 
		if(!version.check()) { 
			System.err.println(String.format("InsertableAlignmentVersion %s failed check.", version.toString()));
			return false; 
		}
		for(InsertableAlignment a : alignments) { 
			if(!a.check()) { 
				System.err.println(String.format("InsertableAlignment %s failed check.", a.toString()));
				return false;
			}
		}
		
		for(InsertableAlignBlock a : blocks) { 
			if(!a.check()) { 
				System.err.println(String.format("InsertableAlignBlock %s failed check.", a.toString()));
				return false;
			}
		}
		System.out.println(String.format("# Alignments: %d", alignments.size()));
		System.out.println(String.format("# Blocks: %d", blocks.size()));
		return version != null && alignments.size() > 0 && blocks.size() > 0;
	}
}
