/*
 * Author: tdanford
 * Date: May 1, 2008
 */
package edu.mit.csail.cgs.utils.datastructures;

import java.util.LinkedList;
import java.io.*;

/**
 * BinarySplayTree
 * @author tdanford
 * 
 * BinarySplayTree is an implementation of binary search trees with the self-adjusting 
 * "splay" maneuver -- these are commonly known as "splay trees."  
 * 
 * This was originally outlined in the paper: 
 * Sleator and Tarjan, "Self-Adjusting Binary Search Trees" (JACM V.32 no.3, July 1985)
 * 
 * All "quoted" comments to the methods in this class are direct quotes taken from that 
 * original paper.  
 *
 * @param <X>
 */
public class BinarySplayTree<X extends Comparable<X>>{

	private X value;
	private BinarySplayTree<X> left, right, parent;
	
	public BinarySplayTree(X v) { 
		value = v;
		parent = left = right = null;
	}
	
	public BinarySplayTree(X v, BinarySplayTree<X> l, BinarySplayTree<X> r) {
		value = v; 
		parent = null; 
		left = l; 
		right = r; 
	}
	
	public BinarySplayTree(X v, BinarySplayTree<X> l, BinarySplayTree<X> r, BinarySplayTree<X> p) { 
		value = v; 
		left = l;
		right = r; 
		parent = p;
	}
	
	public boolean isRoot() { return parent==null; }
	public boolean isLeaf() { return left == null && right == null; }
	public X getValue() { return value; }
	public BinarySplayTree<X> getLeft() { return left; }
	public BinarySplayTree<X> getRight() { return right; }
	
	public boolean isLeftChild() { 
		return !isRoot() && isLeftOf(parent);
	}
	
	public boolean isRightChild() { 
		return !isRoot() && isRightOf(parent);
	}
	
	public boolean isLeftOf(BinarySplayTree<X> pt) { 
		return pt != null && pt==parent && pt.left==this; 
	}
	
	public boolean isRightOf(BinarySplayTree<X> pt) { 
		return pt != null && pt==parent && pt.right==this; 
	}
	
	public void printTree() { 
		printTree(System.out);
	}
	
	public void printTree(PrintStream ps) { 
		printTree(ps, 0);
	}
	
	private void printTree(PrintStream ps, int indent) { 
		ps.println(value.toString());
		if(left != null) { 
			printIndent(ps, indent+1);
			ps.print("-");
			left.printTree(ps, indent+1);
		}
		if(right != null) { 
			printIndent(ps, indent+1);
			ps.print("+");
			right.printTree(ps, indent+1);
		}
	}
	
	private void printIndent(PrintStream ps, int indent) { 
		for(int i = 0; i < indent; i++) { 
			ps.print("  ");
		}		
	}
	
	public int childIndicator() {
		if(isRightChild()) { 
			return 1;
		} else if (isLeftChild()) { 
			return -1; 
		} else { 
			return 0;
		}
	}
	
	public BinarySplayTree<X> search(X v) { 
		BinarySplayTree<X> t = findLast(v);
		return t.value.equals(v) ? t : null;
	}
	
	private BinarySplayTree<X> findLast(X v) { 
		int c = v.compareTo(value);
		if(c == 0) { 
			return this;
		} else if (c < 0) { 
			return left != null ? left.search(v) : this;
		} else { 
			return right != null ? right.search(v) : this;
		}		
	}
	
	private LinkedList<BinarySplayTree<X>> stubPath() { 
		LinkedList<BinarySplayTree<X>> lst = new LinkedList<BinarySplayTree<X>>();
		lst.addLast(this);
		return lst;
	}
	
	public LinkedList<BinarySplayTree<X>> pathingSearch(X v) { 
		int c = v.compareTo(value);
		if(c == 0) { 
			return stubPath();
		} else if (c < 0) { 
			if(left != null) { 
				LinkedList<BinarySplayTree<X>> path = left.pathingSearch(v);
				path.addFirst(this);
				return path;
			} else { 
				return stubPath();
			}
		} else { 
			if(right != null) { 
				LinkedList<BinarySplayTree<X>> path = right.pathingSearch(v);
				path.addFirst(this);
				return path;
			} else { 
				return stubPath();
			}
		}		
	}
	
	public BinarySplayTree<X> findRoot() { 
		return isRoot() ? this : parent.findRoot();
	}
	
	public BinarySplayTree<X> findLargestLeaf() { 
		if(right == null) { 
			return this; 
		} else { 
			return right.findLargestLeaf();
		}
	}

	/**
	 * "To perform access(i, t), we search down from the root of t, looking for i.  
	 * If the search reaches node x containing i, we complete the access by splaying
	 * at x and returning a pointer to x.  If the search reaches a null node, indicating
	 * that i is not in the tree, we complete the access by splaying at the last 
	 * nonnull node reached during the search (the node from which the search ran off
	 * the bottom of the tree) and returning a pointer to null.  If the tree is empty, 
	 * we omit the splaying operation." 
	 * 
	 * @param i
	 * @return
	 */
	public boolean access(X i) { 
		BinarySplayTree<X> leaf = findLast(i);
		leaf.splay();
		return value.equals(i);
	}
	
	/**
	 * "To carry out insert(i, t), we perform split(i, t) and then replace t
	 * by a tree consisting of a new root node containing i, whose left and right
	 * subtrees are the trees t1 and t2 returned by the split."
	 * 
	 * @param i
	 */
	public void insert(X i) {
		BinarySplayTree<X> other = split(i);
		BinarySplayTree<X> newLeft = new BinarySplayTree(value, left, right);
		
		value = i;
		left = newLeft; newLeft.parent = this;
		right = other; right.parent = this;
	}
	
	/**
	 * "To carry out delete(i, t), we perform access(i, t) and then replace t by
	 * the join of its left and right subtrees."
	 * 
	 * @param i
	 */
	public BinarySplayTree<X> delete(X i) { 
		access(i);
		left.parent = right.parent = null;
		left.join(right);
		return left;
	}
	
	/**
	 * "To carry out join(t1, t2), we begin by accessing the largest item, say
	 * i, in t1.  After the access, the root of t1 contains i and thus has a 
	 * null right child.  We complete the join by making t2 the right subtree
	 * of this root and returning the resulting tree." 
	 * 
	 * "In both join and split, we must deal specially with the case of an empty 
	 * input tree (or trees)."
	 * 
	 * @param t2
	 */
	public void join(BinarySplayTree<X> t2) {
		if(!isRoot()) { throw new IllegalArgumentException(); }
		if(!t2.isRoot()) { throw new IllegalArgumentException(); } 
		
		BinarySplayTree<X> t1Largest = findLargestLeaf();
		BinarySplayTree<X> newRoot = t1Largest.splay();
		
		if(newRoot.right != null) { throw new IllegalStateException(); }
		
		newRoot.right = t2;
		t2.parent = newRoot;
	}
	
	/**
	 * "To carry out split(i, t) we perform acces(i, t) and then return the 
	 * two trees formed by breaking either the left link or the right link from 
	 * the new root of t, depending on whether the root contains an item greater 
	 * than i or not greater than i."  
	 * 
	 * "In both join and split, we must deal specially with the case of an empty 
	 * input tree (or trees)."
	 * 
	 * @param i
	 * @return
	 */
	public BinarySplayTree<X> split(X i) {
		if(!isRoot()) { throw new IllegalArgumentException("Can't split() on a non-root node."); }
		
		BinarySplayTree<X> last = findLast(i);
		last.splay();
		
		BinarySplayTree<X> other = null;
		
		if(value.compareTo(i) == 1) { 
			// break left link
			other = left;
			left = null;
			other.parent = null;
		} else { 
			// break right link
			other = right;
			right = null;
			other.parent = null;
		}
		
		return other;
	}
	
	/**
	 * "Splaying step: 
	 * Case 1 (zig): If p(x), the parent of x, is the tree root, rotate the edge
	 * joining x with p(x).  (This case is terminal.)
	 * 
	 * Case 2 (zig-zig): If p(x) is not the root and x and p(x) are both left or
	 * right children, rotate the edge joining p(x) with its grandparent g(x) and 
	 * then rotate the edge joining x with p(x).
	 * 
	 * Case 3 (zig-zag): If p(x) is not the root and x is the left child and p(x) is
	 * the right child (or vice-versa), rotate the edge joining x with p(x), and then 
	 * rotate the edge joining x with the new p(x)." 
	 * 
	 * @return
	 */
	public BinarySplayTree<X> splay() { 
		if(!isRoot()) { 
			if(parent.isRoot()) { 
				rotate();
				return parent;
			} else if (childIndicator() == parent.childIndicator()) {   
				parent.rotate();
				rotate();
				return parent.parent.splay();
			} else { 
				rotate();
				rotate();
				return parent.parent.splay();
			}
		}
		return this;
	}
	
	public void rotate() { 
		if(isLeftChild()) { 
			rotateRight(); 
		} else if (isRightChild()) { 
			rotateLeft();
		}
	}
	
	public void rotateRight() { 
		if(isRoot()) { throw new IllegalArgumentException("Can't rotate right on root."); }
		if(!isLeftOf(parent)) { throw new IllegalArgumentException("Right child can't rotate right."); }

		X temp = value;
		
		parent.left = left; left.parent = parent;
		
		value = parent.value;
		left = right;
		right = parent.right; right.parent = this;
		
		parent.value = temp;
		parent.right = this; 
	}

	public void rotateLeft() { 
		if(isRoot()) { throw new IllegalArgumentException("Can't rotate left on root."); }
		if(!isRightOf(parent)) { throw new IllegalArgumentException("Left child can't rotate left."); }

		X temp = value;
		
		parent.right = right; right.parent = parent;
		
		value = parent.value;
		right = left;
		left = parent.left; left.parent = this;
		
		parent.value = temp;
		parent.left= this; 
	}
}
