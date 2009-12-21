package edu.mit.csail.cgs.clustering.hierarchical;

import java.util.HashSet;
import java.util.Set;

import edu.mit.csail.cgs.clustering.Cluster;

/**
 * @author Timothy Danford
 *
 */
public class ClusterNode<X> implements Cluster<X> {
	
    private Cluster left, right;
    private double weight;
	
    public ClusterNode(Cluster l, Cluster r) { 
        left = l; right = r;
        weight = 1.0;
    }
	
    public ClusterNode(double w, Cluster l, Cluster r) { 
        left = l; right = r;
        weight = w;
    }
	
    public double getWeight() { return weight; }
    public Cluster getLeft() { return left; }
    public Cluster getRight() { return right; }
	
    public Cluster getLeftmost() { 
        if(left instanceof ClusterNode) { 
            return ((ClusterNode)left).getLeftmost();
        } else { 
            return left; 
        }
    }
	
    public Cluster getRightmost() { 
        if(right instanceof ClusterNode) { 
            return ((ClusterNode)right).getRightmost();
        } else { 
            return right;
        }
    }
	
    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.clustering.Cluster#getElements()
     */
    public Set<X> getElements() {
        HashSet<X> ce = new HashSet<X>();
        ce.addAll(left.getElements());
        ce.addAll(right.getElements());
        return ce;
    }

    public int hashCode() { 
        int code = 17;
        code += left.hashCode(); code *= 37;
        code += right.hashCode(); code *= 37;
        return code;
    }
	
    public boolean equals(Object o) { 
        if(!(o instanceof ClusterNode)) { return false; }
        ClusterNode cn = (ClusterNode)o;
        if(!left.equals(cn.left)) { return false; }
        if(!right.equals(cn.right)) { return false; }
        return true;
    }
	
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        sb.append(left.toString());
        sb.append(",");
        sb.append(right.toString());
        sb.append("]");
        return sb.toString();
    }

    public int size() {
        return left.size() + right.size();
    }
	
    public int depth() {
        int ld = 1; 
        if(left instanceof ClusterNode) { ld = ((ClusterNode)left).depth(); }
        int rd = 1; 
        if(right instanceof ClusterNode) { rd = ((ClusterNode)right).depth(); }
        return Math.max(ld, rd) + 1;
    }
}
