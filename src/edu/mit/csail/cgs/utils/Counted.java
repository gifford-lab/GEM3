package edu.mit.csail.cgs.utils;

import java.lang.*;
import java.io.*;
import java.util.*;

public class Counted { 

	private Map fCountMap;

	public Counted() { 
		fCountMap = new HashMap();
	}

	public Counted(Counted c) { 
		fCountMap = new HashMap();
		Iterator itr = c.fCountMap.entrySet().iterator();
		while(itr.hasNext()) { 
			Map.Entry e = (Map.Entry)itr.next();
			fCountMap.put(e.getKey(), new ObjectCount((ObjectCount)e.getValue()));
		}
	}

	public int[] getCountDistribution() { 
		return getCountDistribution(getMinCount(), getMaxCount());
	}

	public int[] getCountDistribution(int start, int end) { 
		if(end < start) { throw new IllegalArgumentException(); }
		int[] counts = new int[end-start+1];
		for(int i = 0; i < counts.length; i++) { 
			counts[i] = 0;
		}

		Iterator itr = fCountMap.values().iterator();
		while(itr.hasNext()) {
			ObjectCount oc = (ObjectCount)itr.next();
			int index = oc.getCount() - start;
			if(index >= 0 && index < counts.length) { counts[index] += 1; }
		}
		return counts;
	}

	public int getMaxCount() { 
		Iterator itr = fCountMap.entrySet().iterator();
		int maxCount = 0;
		while(itr.hasNext()) { 
			Map.Entry e = (Map.Entry)itr.next();
			ObjectCount oc = (ObjectCount)e.getValue();
			if(oc.getCount() > maxCount) { maxCount = oc.getCount(); }
		}
		return maxCount;
	}

	public int getMinCount() { 
		Iterator itr = fCountMap.entrySet().iterator();
		int minCount = -1;
		while(itr.hasNext()) { 
			Map.Entry e = (Map.Entry)itr.next();
			ObjectCount oc = (ObjectCount)e.getValue();
			if(minCount == -1 || oc.getCount() < minCount) { minCount = oc.getCount(); }
		}
		return minCount;
	}

	public void clear() { fCountMap = new HashMap(); }

	public int size() { return fCountMap.size(); }

	public int totalCount() { 
		int total = 0;
		Iterator itr = fCountMap.values().iterator();
		while(itr.hasNext()) { 
			ObjectCount oc = (ObjectCount)itr.next();
			total += oc.getCount();
		}
		return total;
	}

	public void addTo(Counted c) { 
		Iterator itr = fCountMap.values().iterator();
		while(itr.hasNext()) { 
			ObjectCount oc = (ObjectCount)itr.next();
			c.changeCount(oc.getToken(), oc.getCount());
		}
	}

	public Set getTokens() { return fCountMap.keySet(); }
	public boolean hasToken(String tok) { return fCountMap.containsKey(tok); }

	public int getCount(String tok) { 
		if(!fCountMap.containsKey(tok)) { return 0; }
		ObjectCount oc = (ObjectCount)fCountMap.get(tok);
		return oc.getCount();
	}

	public int getDiff(Counted c, String token) { 
		int c1 = getCount(token);
		int c2 = c.getCount(token);
		return c1 - c2;
	}

	public void printDiffs(Counted c) { printDiffs(c, System.out); }
	public void printDiffs(Counted c, PrintStream ps) { 
		Set tokenSet = new HashSet();
		SortedSet diffSet = new TreeSet();
		tokenSet.addAll(getTokens()); tokenSet.addAll(c.getTokens());
		Iterator itr = tokenSet.iterator();
		while(itr.hasNext()) { 
			String tok = (String)itr.next();
			int diff = getDiff(c, tok);
			diffSet.add(new ObjectCount(tok, diff));
		}

		itr = diffSet.iterator();
		while(itr.hasNext()) { 
			ps.println(itr.next());
		}
	}

	public void changeCount(String tok, int count) { 
		if(!hasToken(tok) && count < 1) { 
			throw new IllegalArgumentException();
		}

		if(!hasToken(tok)) { 
			ObjectCount oc = new ObjectCount(tok, count);
			fCountMap.put(tok, oc);
		} else { 
			ObjectCount oc = (ObjectCount)fCountMap.get(tok);
			if(oc.getCount() + count < 0) { 
				throw new IllegalArgumentException();
			}

			oc.changeCount(count);
			if(oc.getCount() == 0) { 
				fCountMap.remove(tok);
			}
		}
	}

	public void printCounts() { printCounts(System.out); }
	public void printCounts(PrintStream ps) { 
		Iterator itr = fCountMap.values().iterator();
		while(itr.hasNext()) { 
			ObjectCount oc = (ObjectCount)itr.next();
			ps.println(oc);	
		}
	}
}

