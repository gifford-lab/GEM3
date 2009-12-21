package edu.mit.csail.cgs.tools.chippet;

import java.io.PrintStream;
import java.util.*;

import edu.mit.csail.cgs.utils.Interval;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class IntervalTree<X> {
    
    private IntervalTree<X> parent;
	private LinkedList<Interval<X>> intervals;
	public IntervalTree<X> left, right;
	private int pivot;
	private Interval totalInterval;
    private long size;
	
	public IntervalTree(Interval intv, Collection<Interval<X>> ints) {
        parent = null;
		intervals = new LinkedList<Interval<X>>(ints);
		left = right = null;
		pivot = 0;
		totalInterval = intv;
        size = ints.size();
	}
	
	public IntervalTree(Interval intv) {
        parent = null;
		totalInterval = intv;
		intervals = new LinkedList<Interval<X>>();
		left = right = null;
		pivot = 0;
        size = 0;
	}
	
	protected IntervalTree(IntervalTree<X> p, Interval intv, Interval<X>[] array, int i1, int i2) {
        parent = p;
		intervals = new LinkedList<Interval<X>>();
		for(int i = i1; i < i2; i++) { intervals.add(array[i]); }
        size = Math.max(0, i2-i1);
		pivot = 0;
		totalInterval = intv;
	}
    
    public IntervalTree<X> getParent() { return parent; }
    
    public void printTreeSummary(PrintStream ps) { printTreeSummary(0, ps); } 
    
    public void printTreeSummary(int indent, PrintStream ps) {
        String istr = "  ";
        String indentstr = "";
        for(int i = 0; i < indent; i++) { indentstr += istr; }
        
        if(intervals != null) { 
            ps.println(indentstr + "* " + totalInterval + " (" + intervals.size() + ")");
			/*
			if(intervals.size() < 10) { 
				for(Interval intv : intervals) { 
					ps.println(indentstr + "  " + intv);
				}
			}
			*/
        } else { 
            ps.println(indentstr + "* " + totalInterval + " (" + size + ")");
            left.printTreeSummary(indent+1, ps);
            right.printTreeSummary(indent+1, ps);
        }
    }
	
	public void addInterval(Interval<X> intv) { 
		if(intervals == null) { 
			if(intv.end < pivot) { left.addInterval(intv); }
			if(intv.start > pivot) { right.addInterval(intv); }
			if(intv.contains(pivot)) { 
				left.addInterval(intv);
				right.addInterval(intv);
			}
		} else { 
			intervals.addLast(intv);
		}
        size += 1;
	}
    
    public long getSize() { return size; }
    
    public long removeValueWithData(X d) { 
        boolean foundValue = false;

        if(intervals != null) {
            Iterator<Interval<X>> itr = intervals.iterator();
            long r = 0;
            while(itr.hasNext()) { 
                Interval<X> intv = itr.next();
                if(intv.data.equals(d)) { 
                    foundValue = true;
                    itr.remove();
                    r += 1;
                }
            }
            
            size -= r;
            return r;
        } else { 
            long r1 = left.removeValueWithData(d);
            long r2 = right.removeValueWithData(d);
            
            size -= (r1 + r2);
            return r1 + r2;
        }
    }

	public void collectLeafIntervals(int pt, LinkedList<Interval<X>> ints) { 
		if(intervals != null) { 
			for(Interval<X> intv : intervals) { 
				if(intv.contains(pt)) { 
					ints.addLast(intv);
				}
			}
		} else { 
			if(pt < pivot) { 
				left.collectLeafIntervals(pt, ints);
			} else if (pt > pivot) { 
				right.collectLeafIntervals(pt, ints);
			} else { 
				left.collectLeafIntervals(pt, ints);
				right.collectLeafIntervals(pt, ints);
			}
		}
	}
	
	public void collectLeafIntervals(Interval query, LinkedList<Interval<X>> ints) { 
		if(intervals != null) { 
			for(Interval<X> intv : intervals) { 
				if(intv.overlaps(query)) { 
					ints.addLast(intv);
				}
			}
		} else { 
			if(query.end < pivot) { 
				left.collectLeafIntervals(query, ints);
			} else if (query.start > pivot) { 
				right.collectLeafIntervals(query, ints);
			} else { 
				Pair<Interval,Interval> sp = query.split(pivot);
				left.collectLeafIntervals(sp.getFirst(), ints);
				right.collectLeafIntervals(sp.getLast(), ints);
			}
		}		
	}

	public void collectLeafIntervals(LinkedList<Interval<X>> ints) { 
		if(intervals != null) { 
			ints.addAll(intervals);
		} else { 
			left.collectLeafIntervals(ints);
			right.collectLeafIntervals(ints);
		}
	}
	
	public boolean isLeaf() { 
		return intervals != null;
	}
	
    public void splitForDepth(int depth) { 
        splitForDepth(depth, 1000);
    }
    
	public void splitForDepth(int depth, long maxSize) {
	    if(size > maxSize) { 
	        if(!isLeaf()) {
	            left.splitForDepth(depth, maxSize);
	            right.splitForDepth(depth, maxSize);
	        } else { 
	            findPivot();

	            if(depth > 1) { 
	                left.splitForDepth(depth-1, maxSize);
	                right.splitForDepth(depth-1, maxSize);
	            }
	        }
	    }
	}
    
    private void findPivot() {  
        Interval<X>[] array = new Interval[intervals.size()];
        int i = 0;
        for(Interval<X> intv : intervals) { array[i++] = intv; }
        Arrays.sort(array);
        
        PivotData d = findPivotFromArray(array);        
        splitOnPivot(d, array);
    }
	
    private void splitOnPivot(PivotData d, Interval<X>[] array) { 
        pivot = d.pivotValue;
        intervals = null;
        left = new IntervalTree<X>(this, new Interval(totalInterval.start, pivot), array, 0, d.rightIndex);
        right = new IntervalTree<X>(this, new Interval(pivot, totalInterval.end), array, d.leftIndex+1, array.length);        
    }
    
	private PivotData findPivotFromArray(Interval[] array) {
		if(array == null || array.length < 1) { 
			return new PivotData(totalInterval.getMidpoint(), 0, 0);
		}

		Arrays.sort(array, new StartComparator());
		int rightIndex = array.length/2;
		
		Arrays.sort(array, 0, rightIndex, new EndComparator());
		int leftIndex = rightIndex;
		
		int pivotValue = array[rightIndex].start;

		while(leftIndex >= 0 && array[leftIndex].end >= pivotValue) { 
			leftIndex -= 1;
		}
		
		return new PivotData(pivotValue, leftIndex, rightIndex);
	}
	
	public Interval getTotalInterval() { 
		return totalInterval;
	}
	
	public int getPivot() { 
		return pivot;
	}
    
    public Iterator<Interval<X>> getIterator() { return new IntervalTreeIterator<X>(this); }

    private static class IntervalTreeIterator<X> implements Iterator<Interval<X>> {
        
        private Iterator<Interval<X>> itr;
        private Stack<IntervalTree<X>> treestack;
        
        public IntervalTreeIterator(IntervalTree<X> c) { 
            itr = null;
            treestack = new Stack<IntervalTree<X>>();
            
            itr = new EmptyIterator<Interval<X>>();
            treestack.push(c);
            findNextIterator();
        }
        
        private void findNextIterator() { 
            while(!itr.hasNext() && !treestack.isEmpty()) { 
                IntervalTree<X> nextit = treestack.pop();
                //System.out.println(nextit.totalInterval);
                if(nextit.intervals != null) { 
                    itr = nextit.intervals.iterator();
                } else { 
                    treestack.push(nextit.right);
                    treestack.push(nextit.left);
                }
            }
        }

        public boolean hasNext() {
            return itr != null && itr.hasNext();
        }

        public Interval<X> next() {
            Interval<X> ret = itr.next();
            if(!itr.hasNext()) { 
                findNextIterator();
            }
            return ret;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        } 
    }    
}

class PivotData { 
	public int pivotValue;
	public int leftIndex, rightIndex;
	
	public PivotData(int v, int l, int r) { 
		pivotValue = v;
		leftIndex = l;
		rightIndex = r;
	}
}

class StartComparator implements Comparator<Interval> { 
	
	public int compare(Interval i1, Interval i2) { 
		if(i1.start < i2.start) { return -1; }
		if(i1.start > i2.start) { return 1; }
		if(i1.end < i2.end) { return -1; }
		if(i1.end > i2.end) { return 1; }
		return 0;
	}
}

class EndComparator implements Comparator<Interval> {
	public int compare(Interval i1, Interval i2) {
		if(i1.end < i2.end) { return -1; }
		if(i1.end > i2.end) { return 1; }
		if(i1.start < i2.start) { return -1; }
		if(i1.start > i2.start) { return 1; }
		return 0;
	}
	
}
