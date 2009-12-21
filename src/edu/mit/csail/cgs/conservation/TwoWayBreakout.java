package edu.mit.csail.cgs.conservation;

import java.util.*;
import edu.mit.csail.cgs.utils.*;

/** 
 * @author tdanford
 *
 * TwoWayBreakout : creates BindingPartition objects from datasets, 
 * for any given ExptDescriptor object.
 */
public class TwoWayBreakout {
	
	private static String fversion;
	
	static { 
		fversion = "GO:5_23_2006";
	}
	
	private char s1, s2;
	private ConservationDataset d1, d2;
	private BindingOptions opts;
	private GeneMap m1, m2;
	private SetTools<String> tools;
	
	public TwoWayBreakout(BindingOptions opts,
			char s1, ConservationDataset d1, GeneMap m1,  
			char s2, ConservationDataset d2, GeneMap m2) { 

		this.s1 = s1; this.s2 = s2;
		this.opts = opts;
		this.d1 = d1; this.d2 = d2;
		this.m1 = m1; this.m2 = m2;
		
		tools = new SetTools<String>();
	}
	
	public TwoWayBreakout getReverseBreakout() { 
		return new TwoWayBreakout(opts, s2, d2, m2, s1, d1, m1);
	}
	
	public BindingPartition createPartition(String shortTitle, 
			ExptDescriptor ed1, ExptDescriptor ed2) {
		
		BindingPartition part = new BindingPartition(shortTitle, ed1.getName(), fversion); 
		
		Set<String> total = d1.getIDs();
		HashSet<String> bound = new HashSet<String>();
		HashSet<String> shadow = new HashSet<String>();
		
		for(String id : total) { 
			if(d1.isBound(id, ed1, opts)) { bound.add(id); }
			
			Set<String> mapped = m1.mapID(id);
			Iterator<String> itr = mapped.iterator();
			boolean foundMatch = false;
			while(itr.hasNext() && !foundMatch) { 
				String mid = itr.next();
				if(d2.isBound(mid, ed2, opts)) { 
					foundMatch = true;
				}
			}
			
			if(foundMatch) { shadow.add(id); }
		}
		
		String consKey = createKey(true, true);
		String boundKey = createKey(true, false);
		String shadKey = createKey(false, true);
		String noneKey = createKey(false, false);
		
		part.addBlock(consKey, tools.intersection(bound, shadow));
		part.addBlock(boundKey, tools.subtract(bound, shadow));
		part.addBlock(shadKey, tools.subtract(shadow, bound));
		part.addBlock(noneKey, tools.subtract(total, tools.union(bound, shadow)));
		
		return part;
	}
	
	public String createKey(boolean b1, boolean b2) { 
		String st1 = String.valueOf(b1 ? Character.toUpperCase(s1) : Character.toLowerCase(s1));
		String st2 = String.valueOf(b2 ? Character.toUpperCase(s2) : Character.toLowerCase(s2));
		return st1 + st2;
	}
}
