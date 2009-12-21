package edu.mit.csail.cgs.ewok.verbs;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.SetTools;

public class UniqueishGeneFilter<X extends Gene> implements Filter<X, X> {

	private int startBuffer;
	private int endBuffer;
	private Map<Pair<String,Integer>,Set<X>> startSeen;
	private Map<Pair<String,Integer>,Set<X>> endSeen;
	private SetTools<X> st;
	
	public UniqueishGeneFilter(int startBuffer, int endBuffer) {
		this.startBuffer = startBuffer;
		this.endBuffer = endBuffer;
		startSeen = new HashMap<Pair<String,Integer>,Set<X>>();
		endSeen = new HashMap<Pair<String,Integer>,Set<X>>();
		this.st = new SetTools<X>();
	}

	public X execute(X a) {
		int start, end;
		switch(a.getStrand()) { 
		default:
		case '+':
			start = a.getStart();
			end = a.getEnd();
		case '-':
			start = a.getEnd();
			end = a.getStart();
		}
		Set<X> startSet = new HashSet<X>();
		if (startBuffer!=-1) {
			for (int l=start-startBuffer; l<=start+startBuffer; l++) {
				Pair<String,Integer> tmpPair = new Pair<String,Integer>(a.getChrom(),l);
				if (startSeen.containsKey(tmpPair)) {
					startSet.addAll(startSeen.get(tmpPair));
				}
				if (l==3621296) {
					if (startSeen.containsKey(tmpPair)) {
						//System.err.println(startSeen.get(tmpList).size());
					} else {
						//System.err.println("empty");
					}
				}
			}
		}
		Set<X> endSet = new HashSet<X>();
		if (endBuffer!=-1) {
			for (int r=end-endBuffer; r<=end+endBuffer; r++) {
				Pair<String,Integer> tmpPair = new Pair<String,Integer>(a.getChrom(),r);
				if (endSeen.containsKey(tmpPair)) {
					endSet.addAll(endSeen.get(tmpPair));
				}
			}
		}
		if ((startBuffer==-1 && endSet.size()>0) ||
				(endBuffer==-1 && startSet.size()>0) ||
				st.intersects(startSet,endSet)) {
			//System.err.println(a);
			return null;
		}
		Pair<String,Integer> tmpPair = new Pair<String,Integer>(a.getChrom(),start);
		if (!startSeen.containsKey(tmpPair)) {
			startSeen.put(tmpPair, new HashSet<X>());
		}
		startSeen.get(tmpPair).add(a);
		tmpPair = new Pair<String,Integer>(a.getChrom(),end);
		if (!endSeen.containsKey(tmpPair)) {
			endSeen.put(tmpPair, new HashSet<X>());
		}
		endSeen.get(tmpPair).add(a);
		return a;
	}

}
