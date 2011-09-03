package edu.mit.csail.cgs.utils;

import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.datastructures.Taxonomy;

public class SetTools<X> {
	
	public SetTools() { 
		
	}
	
	public boolean isEqual(Set<X> first, Set<X> second) { 
		if(first.size() != second.size()) { return false; }
		for(X v : first) { if(!second.contains(v)) { return false; } }
		return true;
	}
	
	public Set<X> union(Set<X> first, Set<X> second) { 
		HashSet<X> total = new HashSet<X>(first);
		total.addAll(second);
		return total;
	}
    
    public Set<X> intersection(Set<X> first, Set<X> second) { 
        Set<X> inter = new HashSet<X>();
        Set<X> iter, other;
        if(first == null || second == null) { throw new IllegalArgumentException("Both the first and second sets should not equal to null."); }
        if(first.size() <= second.size()) { iter = first;  other = second; } 
        else 							  { iter = second; other = first;  }
        for(X v : iter) { 
            if(other.contains(v)) { inter.add(v); }
        }
        return inter;
    }
	
	public boolean intersects(Set<X> first, Set<X> second) { 
		Set<X> test = first, other = second;
		if(second.size() < first.size()) { test = second; other = first; }
		for(X v : test) { if(other.contains(v)) { return true; } }
		return false; 
	}

	public boolean contains(Set<X> outer, Set<X> inner) { 
		for(X v : inner) { if(!outer.contains(v)) { return false; } }
		return true; 
	}
    
	/**
	 * Get the elements in first set but not in second set
	 * @param first
	 * @param second
	 * @return
	 */
    public Set<X> subtract(Set<X> first, Set<X> second) { 
        HashSet<X> s = new HashSet<X>();
        for(X v : first) { if(!second.contains(v)) { s.add(v); } }
        return s;
    }
    
    public Collection<Pair<LinkedList<String>,X>> collectAddressedItems(Taxonomy<X> tax) { 
        LinkedList<Pair<LinkedList<String>,X>> lst = 
            new LinkedList<Pair<LinkedList<String>,X>>();
        
        for(X v : tax.getImmediateElements()) { 
            LinkedList<String> empty = new LinkedList<String>();
            Pair<LinkedList<String>,X> p = new Pair<LinkedList<String>,X>(empty, v);
            lst.addLast(p);
        }
        
        for(String addr : tax.getAddrs()) { 
            Taxonomy<X> subTax = tax.getSubTaxonomy(addr);
            Collection<Pair<LinkedList<String>,X>> subRes = collectAddressedItems(subTax);
            for(Pair<LinkedList<String>,X> p : subRes) { 
                p.getFirst().addLast(addr);
                lst.addLast(p);
            }
        }
            
        return lst;
    }
}
