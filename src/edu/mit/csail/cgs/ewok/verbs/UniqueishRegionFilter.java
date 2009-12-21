package edu.mit.csail.cgs.ewok.verbs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.utils.SetTools;

public class UniqueishRegionFilter<X extends Region> implements Filter<X, X> {

	private int buffer;
	private Map<List<String>,Set<X>> leftSeen;
	private Map<List<String>,Set<X>> rightSeen;
	private SetTools<X> st;
	
	public UniqueishRegionFilter(int buffer) {
		this.buffer = buffer;
		this.leftSeen = new HashMap<List<String>,Set<X>>();
		this.rightSeen = new HashMap<List<String>,Set<X>>();
		this.st = new SetTools<X>();
	}
	
	public X execute(X a) {
		List<String> tmpList = new ArrayList<String>();
		tmpList.add(a.getChrom());
		Set<X> leftSet = new HashSet<X>();
		for (int l=a.getStart()-buffer; l<=a.getStart()+buffer; l++) {
			tmpList.add((new Integer(l)).toString());
			if (leftSeen.containsKey(tmpList)) {
				leftSet.addAll(leftSeen.get(tmpList));
			}
			tmpList.remove(1);
		}
		Set<X> rightSet = new HashSet<X>();
		for (int r=a.getEnd()-buffer; r<=a.getEnd()+buffer; r++) {
			tmpList.add((new Integer(r)).toString());
			if (rightSeen.containsKey(tmpList)) {
				rightSet.addAll(rightSeen.get(tmpList));
			}
			tmpList.remove(1);
		}
		if (st.intersects(leftSet,rightSet)) {
			return null;
		}
		tmpList.add((new Integer(a.getStart())).toString());
		if (!leftSeen.containsKey(tmpList)) {
			leftSeen.put(tmpList, new HashSet<X>());
		}
		leftSeen.get(tmpList).add(a);
		tmpList.remove(1);
		tmpList.add((new Integer(a.getEnd())).toString());
		if (!rightSeen.containsKey(tmpList)) {
			rightSeen.put(tmpList, new HashSet<X>());
		}
		rightSeen.get(tmpList).add(a);
		return a;
	}

}
