package edu.mit.csail.cgs.utils;


/**
 * @author tdanford
 */
public class Interval<X> implements Comparable<Interval> {

	public int start, end; // inclusive coordinates.
    public X data;
	
	public Interval(int s, int e) { start = s; end = e; data = null; }
    public Interval(int s, int e, X d) { start = s; end = e; data = d; }
	
	public int hashCode() { 
		int code = 17;
		code += start; code *= 37;
		code += end; code *= 37;
        if(data != null) { code += data.hashCode(); code *= 37; }
		return code;
	}
	
	public String toString() { return "[" + start + "," + end + "]"; }
    
    public Interval<X> expand(int w) { 
        return new Interval<X>(start-w, end+w, data);
    }
    
    public int getWidth() { return end-start+1; }
    public int getMidpoint() { return (end+start)/2; }
	
	public Pair<Interval<X>,Interval<X>> split(int pt) { 
		if(!contains(pt)) { throw new IllegalArgumentException(); }
		return new Pair<Interval<X>,Interval<X>>(new Interval<X>(start, pt, data), new Interval<X>(pt, end, data));
	}
	
	public Interval<X> intersection(Interval intv) { 
		if(!overlaps(intv)) { throw new IllegalArgumentException(); }
		return new Interval<X>(Math.max(start, intv.start), Math.min(end, intv.end), data);
	}
	
	public boolean equals(Object o) { 
		if(!(o instanceof Interval)) { return false; }
		Interval i = (Interval)o;
        
        if(data != null || i.data != null) { 
            if(data != null && i.data != null && !data.equals(i.data)) { 
                return false; 
            }
            return false;
        }
		return start == i.start && end == i.end;
	}
	
	public boolean contains(int pt) { return start <= pt && end >= pt; }
	
	public boolean contains(Interval i) { 
		return start <= i.start && end >= i.end;
	}
	
	public boolean overlaps(Interval i) { 
		return contains(i.start) || i.contains(start);
	}
	
	public int overlap(Interval i) {
		if(!overlaps(i)) { 
			return 0; 
		} else if(contains(i)) { 
			return i.getWidth();
		} else if(i.start < start) {   
			return i.end-start+1;
		} else { 
			return end-i.start+1;
		}
	}
	
	public int distance(Interval i) { 
		if(overlaps(i)) { 
			return 0;
		} else { 
			if(i.end < start) { 
				return start-i.end;
			} else { 
				return i.start-end;
			}
		}
	}

	public int compareTo(Interval i) {
		if(start < i.start) { return -1; }
		if(start > i.start) { return 1; }
		if(end < i.end) { return -1; }
		if(end > i.end) { return 1; }
		return 0;
	}
	
}
