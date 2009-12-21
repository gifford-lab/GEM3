package edu.mit.csail.cgs.utils.datastructures;

import java.io.*;
import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.utils.EventSource.Default;

public class TaxonomyImpl<X> implements Taxonomy<X>, Saveable, Listener<TaxonomyEvent> {
	
	private Set<X> elmts;
	private Map<String,TaxonomyImpl<X>> subTaxa;
	private EventSource.Default<TaxonomyEvent> src;
	
	public TaxonomyImpl(DataInputStream dis, Factory<X> loadingFactory) throws IOException { 
		int numElmts = dis.readInt();
		int numSubTaxa = dis.readInt();
		src = new EventSource.Default<TaxonomyEvent>();
		elmts = new HashSet<X>();
		subTaxa = new HashMap<String,TaxonomyImpl<X>>();
		
		for(int i = 0; i < numElmts; i++) { 
			elmts.add(loadingFactory.createObject());
		}
		
		for(int i = 0; i < numSubTaxa; i++) { 
			String k = dis.readUTF();
			TaxonomyImpl<X> ti = new TaxonomyImpl<X>(dis, loadingFactory);
			ti.addEventListener(this);
			subTaxa.put(k, ti);
		}
	}
	
	public TaxonomyImpl() { 
		elmts = new HashSet<X>();
		subTaxa = new HashMap<String,TaxonomyImpl<X>>();
		src = new EventSource.Default<TaxonomyEvent>();
	}
	
	public TaxonomyImpl(Taxonomy<X> t) { 
		src = new EventSource.Default<TaxonomyEvent>();
		if(t != null) { 
			elmts = new HashSet<X>(t.getImmediateElements());
			subTaxa = new HashMap<String,TaxonomyImpl<X>>();
			for(String k : t.getAddrs()) { 
				subTaxa.put(k, new TaxonomyImpl<X>(t.getSubTaxonomy(k)));
				subTaxa.get(k).addEventListener(this);
			}
		} else { 
			elmts = new HashSet<X>();
			subTaxa = new HashMap<String,TaxonomyImpl<X>>();
		}
	}
    
	public void save(DataOutputStream dos) throws IOException { 
		dos.writeInt(elmts.size()); dos.writeInt(subTaxa.size());
		for(X v : elmts) { 
			Saveable s = (Saveable)v;
			s.save(dos);
		}
		for(String k : subTaxa.keySet()) {
			dos.writeUTF(k);
			subTaxa.get(k).save(dos);
		}
	}
	
	public void recursiveRemoveElement(X v) { 
		if(elmts.contains(v)) { 
			elmts.remove(v); 
			fireElementRemoved(v);
		}
		for(String k : subTaxa.keySet()) { 
			subTaxa.get(k).recursiveRemoveElement(v);
		}
	}
	
	public void fireElementRemoved(Object elmt) { 
		//System.out.println("Element Removed: " + elmt);
		TaxonomyEvent te = new TaxonomyEvent(this, TaxonomyEvent.ELEMENT_REMOVED, elmt);
		src.fireEvent(te);		
	}
	
	public void fireElementAdded(Object elmt) { 
		//System.out.println("Element Added: " + elmt);
		TaxonomyEvent te = new TaxonomyEvent(this, TaxonomyEvent.ELEMENT_ADDED, elmt);
		src.fireEvent(te);		
	}
	
	public void fireTaxonRemoved(String addr) { 
		//System.out.println("Taxon Removed: " + addr);
		TaxonomyEvent te = new TaxonomyEvent(this, TaxonomyEvent.TAXON_REMOVED, addr);
		te.addAddr(addr);
		src.fireEvent(te);		
	}
	
	public void fireTaxonAdded(String addr) { 
		//System.out.println("Taxon Added: " + addr);
		TaxonomyEvent te = new TaxonomyEvent(this, TaxonomyEvent.TAXON_ADDED, addr);
		te.addAddr(addr);
		src.fireEvent(te);		
	}
	
	public void fireElementChanged(Object elmt) { 
		//System.out.println("Element Changed: " + elmt);
		if(!elmts.contains(elmt)) { throw new IllegalArgumentException(); }
		TaxonomyEvent te = new TaxonomyEvent(this, TaxonomyEvent.ELEMENT_CHANGED, elmt);
		src.fireEvent(te);
	}
	
	public void eventRegistered(TaxonomyEvent te) { 
		Object source = te.getSource();
		String addr = null;
		Iterator<String> kitr = subTaxa.keySet().iterator();
		while(kitr.hasNext() && addr == null) {
			String k = kitr.next();
			if(subTaxa.get(k) == source) { addr = k; }
		}
		if(addr != null) { 
			TaxonomyEvent te2 = new TaxonomyEvent(this, te, addr); 
			src.fireEvent(te2);
			//System.out.println("** Passing: " + te2.toString());
		} else { 
			//System.out.println("## couldn't find addr.");
			src.fireEvent(te);
		}
	}
	
	public void print(PrintStream ps) { 
		indentedPrint(ps, 0);
	}

	private void indentedPrint(PrintStream ps, int indents) { 
		for(X v : elmts) { 
			for(int i = 0; i < indents; i++) { 
				ps.print("\t");
			}
			ps.println("- " + v.toString());
		}
		for(String k : subTaxa.keySet()) { 
			for(int i = 0; i < indents; i++) { 
				ps.print("\t");
			}
			ps.println("* " + k);
			subTaxa.get(k).indentedPrint(ps, indents+1);
		}
	}
	
	public int size() { 
		int c = elmts.size();
		for(String k : subTaxa.keySet()) { 
			c += subTaxa.get(k).size();
		}
		return c;
	}
	
	public int getNumElements() { return elmts.size(); }
	public int getNumSubTaxonomies() { return subTaxa.size(); }
	
	public void addElement(X k) { 
		if(!elmts.contains(k)) { 
			elmts.add(k);
			fireElementAdded(k);
		}
	}
	
	public void addElement(Collection<String> addrs, X e) {
		if(addrs.size() == 0) { 
			addElement(e);
		} else { 
			LinkedList<String> rest = new LinkedList<String>(addrs);
			String k = rest.removeFirst();
			if(!subTaxa.containsKey(k)) { 
				subTaxa.put(k, new TaxonomyImpl<X>());
				subTaxa.get(k).addEventListener(this);
			}
			subTaxa.get(k).addElement(rest, e);
		}
	}
	
	public void addTaxonomy(String addr, Taxonomy<X> tax) {
		if(subTaxa.containsKey(addr)) { 
			throw new IllegalArgumentException("addr already exists: " + addr);
		}
		subTaxa.put(addr, new TaxonomyImpl<X>(tax));
		subTaxa.get(addr).addEventListener(this);
		fireTaxonAdded(addr);
	}
	
	public Collection<X> getImmediateElements() { 
		return elmts;
	}
	
	public Collection<X> getAllElements() { 
		Set<X> total = new HashSet<X>(elmts);
		for(String k : subTaxa.keySet()) { 
			total.addAll(subTaxa.get(k).getAllElements()); 
		}
		return total;
	}
	
	public Taxonomy<X> getSubTaxonomy(Collection<String> addr) {
		if(addr.size() == 0) { return this; }
		LinkedList<String> lst = new LinkedList<String>(addr);
		String f = lst.removeFirst();
		Taxonomy<X> sub = getSubTaxonomy(f);
		if(sub == null) { return null; }
		return sub.getSubTaxonomy(lst);
	}

	public TaxonomyImpl<X> getSubTaxonomy(String addr) { 
		if(subTaxa.containsKey(addr) && subTaxa.get(addr)==null) { 
			throw new IllegalArgumentException(); 
		}
		return subTaxa.get(addr);
	}
	
	public Set<String> getAddrs() { return subTaxa.keySet(); }
	public boolean hasAddr(String s) { return subTaxa.containsKey(s); }
	
	public Object clone() { 
		TaxonomyImpl<X> newT = new TaxonomyImpl<X>(this);
        return newT;
    }
	
	public Taxonomy<X> combine(Taxonomy<X> t) { 
		TaxonomyImpl<X> newT = new TaxonomyImpl<X>();
		
		newT.elmts.addAll(elmts);		
		for(String k : subTaxa.keySet()) { 
			if(t.hasAddr(k)) { 
				TaxonomyImpl<X> newImpl = new TaxonomyImpl<X>(subTaxa.get(k));
				newImpl = new TaxonomyImpl<X>(t.getSubTaxonomy(k).combine(newImpl));
				newT.subTaxa.put(k, newImpl);
				newImpl.addEventListener(newT);
			} else { 
				newT.subTaxa.put(k, new TaxonomyImpl<X>(subTaxa.get(k)));
				newT.subTaxa.get(k).addEventListener(newT);
			}
		}

		newT.elmts.addAll(t.getImmediateElements());
		for(String k : t.getAddrs()) { 
			if(!subTaxa.containsKey(k)) { 
				newT.subTaxa.put(k, new TaxonomyImpl<X>(t.getSubTaxonomy(k)));
				newT.subTaxa.get(k).addEventListener(newT);
			}
		}
		
		return newT;
	}

	public void addEventListener(Listener<TaxonomyEvent> el) {
		src.addEventListener(el);
	}

	public void removeEventListener(Listener<TaxonomyEvent> el) {
		src.removeEventListener(el);
	}

    public boolean hasListeners() {
        return src.hasListeners();
    }
}
