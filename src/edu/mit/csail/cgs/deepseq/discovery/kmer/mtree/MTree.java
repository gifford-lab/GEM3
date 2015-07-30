package edu.mit.csail.cgs.deepseq.discovery.kmer.mtree;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;

import edu.mit.csail.cgs.deepseq.discovery.kmer.KMAC1;
import edu.mit.csail.cgs.deepseq.discovery.kmer.Kmer;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.utils.sequence.SequenceUtils;


public class MTree {
	
	private MTreeNode root;
	private int capacity;
	private int size;
	private ArrayList<Kmer> data;
	private int currentEntry;
	
	class MTreeNode {
		
		/**
		 * This is the class which represent Node objects. A MTreeNode contains a variety of TreeObjects
		 * and can have one TreeObject as a parent, or null as a parent if it represents the root node.
		 * numObjects represents the number of objects that are in the subtree formed with this node as
		 * root, and does NOT simply refer to the size of nro.
		 * 
		 * @param isLeaf true if this node is not the root, false otherwise
		 * @param nro node routing objects, the objects that this node contains
		 * @param parent null if this node is the root, the object which refers to this node otherwise
		 * @param numObjects the number of objects in the entire tree with this node as root
		 */
		
		private boolean isLeaf;
		private ArrayList<TreeObject> nro;
		private TreeObject parent;
		private int numObjects;
		
		public MTreeNode(boolean l, TreeObject p) {
			isLeaf = l;
			parent = p;
			nro = new ArrayList<TreeObject>();
			numObjects = 0;
		}
		
		public MTreeNode(boolean l, TreeObject p, ArrayList<TreeObject> objects) {
			isLeaf = l;
			parent = p;
			nro = objects;
			numObjects = 0;
		}
		
		public void addTO(TreeObject to) {
			if (to.getContainer() != null) {
				MTreeNode container = to.getContainer();
				int containerSize = container.getObjects().size();
				for (int i = 0; i < containerSize; i++) {
					if (container.getObjects().get(i).getData().getIndex() == (to.getData().getIndex())) {
						container.getObjects().remove(i);
						break;
					}
				}
			}
			nro.add(to);
			to.setContainer(this);
			if (this.isLeaf && (to.getChild() != null && to.getChild().getObjects().size() != 0)) {
				this.setLeaf(false);
			}
		}
		
		public boolean isLeaf() {
			return isLeaf;
		}
		
		public void setLeaf(boolean b) {
			isLeaf = b;
		}
		
		public ArrayList<TreeObject> getObjects() {
			return nro;
		}
		
		public TreeObject getParent() {
			if (parent != null) {
				return parent;
			}
			return null;
		}
		
		public void setParent(TreeObject o) {
			parent = o;
		}
		
		public ArrayList<Kmer> getData() {
			ArrayList<Kmer> data = new ArrayList<Kmer>();
			for (TreeObject o: this.nro) {
				data.add(o.getData());
			}
			return data;
		}
		
		public int getNumObjects() { 
			return numObjects;
		}
		
		public void setNumObjects() {
			int total = this.nro.size();
			if (!this.isLeaf()) {
				for (TreeObject o: this.nro) {
					if (o.getChild() != null) {
						total += o.getChild().getNumObjects();
					}
				}
			}
			numObjects = total;
		}
		
		public void recursiveSetNumObjects() {
			this.setNumObjects();
			for (TreeObject o: nro) {
				if (o.getChild() != null && o.getChild().getObjects().size() > 0) {
					o.getChild().recursiveSetNumObjects();
				}
			}
		}
		
		public void clearTOs() {
			nro.clear();
		}
		
		public void recursivePrint() {
			for (TreeObject o: nro) {
				if (o.getChild() == null || o.getChild().getObjects().size() == 0) {
					if (this.getParent() == null) {
						System.out.println(o.getData().getIndex() + " child of null" + " has radius " + o.getR());
					}
					else {
						System.out.println(o.getData().getIndex() + " child of " + this.getParent().getData().getIndex() + " has radius " + o.getR());
					}
				}
				else {
					o.getChild().recursivePrint();
					if (this.getParent() == null) {
						System.out.println(o.getData().getIndex() + " child of null" + " has radius " + o.getR());
					}
					else {
						System.out.println(o.getData().getIndex() + " child of " + this.getParent().getData().getIndex() + " has radius " + o.getR());
					}
				}
			}
		}
		
		public void recursiveRI() {
			for (TreeObject o: this.getObjects()) {
				o.RI();
				if (o.getChild() != null && o.getChild().getObjects().size() > 0) {
					o.getChild().recursiveRI();
				}
			}
		}
		
		public ArrayList<TreeObject> getLeaves() {
			ArrayList<TreeObject> leaves = new ArrayList<TreeObject>();
			for (TreeObject o: this.getObjects()) {
				if (o.getChild() == null || o.getChild().getObjects().size() == 0) {
					leaves.add(o);
				}
				else {
					leaves.addAll(o.getChild().getLeaves());
				}
			}
			return leaves;
		}
		
		public ArrayList<MTreeNode> getLeafNodes() {
			ArrayList<MTreeNode> leaves = new ArrayList<MTreeNode>();
			if (this.isLeaf()) {
				leaves.add(this);
			}
			else {
				for (TreeObject o: this.getObjects()) {
					if (o.getChild() != null && o.getChild().getObjects().size() > 0) {
						leaves.addAll(o.getChild().getLeafNodes());
					}
				}
			}
			return leaves;
		}
	}
	
	class TreeObject {
		
		/**
		 * This is the class that represents one of the objects contained in the tree. It 
		 * doubles as a routing object if it has a child node associated and as a data object which
		 * contains a Kmer. It has a radius r such that every descendant (of any number of generations)
		 * of the node lies within the ball of radius r centered at the data point (Kmer) contained by
		 * the object, measured by ycDistance. 
		 * 
		 * @param container this is the containing MTreeNode which this object belongs to.
		 * @param radius this is the aforementioned covering radius so that the object contains its 
		 * subtree, if any.
		 * @param child the root of the tree which this object covers, if any.
		 * @param data the actual data object stored in this tree, a Kmer.
		 * @param pnRatio used for the specific data set, calculating ycDistance
		 * @param cut used for the specific data set, calculating ycDistance
		 */
		
		private MTreeNode container;
		private double radius; // covering radius
		private MTreeNode child; // covering tree
		private Kmer data;
		
		public TreeObject(MTreeNode co, MTreeNode c, Kmer d, double r) {
			container = co;
			child = c;
			data = d;
			radius = r;
		}
		
		public TreeObject(MTreeNode co, Kmer d) {
			container = co;
			child = null;
			data = d;
			radius = 0;
		}
		
		public Kmer getData() {
			return data;
		}
		
		public void setData(Kmer d) {
			data = d;
		}
		
		public double getR() {
			return radius;
		}
		
		public void setR(double r) {
			radius = r;
		}
		
		public MTreeNode getChild() {
			return child;
		}
		
		public void setChild(MTreeNode c) {
			this.child = c;
			c.setParent(this);
		}
		
		public MTreeNode getContainer() {
			return container;
		}
		
		public void setContainer(MTreeNode newContainer) {
			container = newContainer;
		}
		
		public void RI() {
			if (this.getChild() == null || this.getChild().getObjects().size() == 0) {
				if (this.getR() != 0) {
					System.out.println("help");
				}
				System.out.println(this.getData().getIndex() + " has radius 0, as expected.");
			}
			else {
				if (this.getR() < this.maxR(this.getChild())) {
					System.out.println("help");
				}
				System.out.println(this.getData().getIndex() + " has radius " + this.getR() + " which beats " + this.maxR(this.getChild()) + ", as expected.");
			}
		}
		
		public double maxR(MTreeNode n) {
			double maxR = 0;
			for (TreeObject o: n.getObjects()) {
				maxR = Math.max(maxR, KMAC1.ycDistance(this.getData(), o.getData()));
				if (o.getChild() != null && o.getChild().getObjects().size() > 0) {
					maxR = Math.max(maxR, this.maxR(o.getChild()));
				}
			}
			return maxR;
		}
		
		public void setParentsR() {
			MTreeNode currentContainer = this.getContainer();
			TreeObject currentAncestor = currentContainer.getParent();
			while (currentAncestor != null) {
				currentAncestor.setR(Math.max(currentAncestor.getR(), KMAC1.ycDistance(this.getData(), currentAncestor.getData())));
				currentContainer = currentAncestor.getContainer();
				currentAncestor = currentContainer.getParent();
			}
		}
	}
	
	public MTree(int k, int s) {
		capacity = k;
		root = new MTreeNode(true, null);
		size = s;
		data = new ArrayList<Kmer>();
		for (int i = 0; i < size; i++) {
			data.add(new Kmer());
		}
		currentEntry = -1;
	}
	
	public MTreeNode getRoot() {
		return root;
	}
	
	public void insertNode2(MTreeNode n, TreeObject entry) {
		data.add(entry.getData().getIndex(), entry.getData());
		if (n.isLeaf()) {
			entry.setContainer(n);
			if (n.getObjects().size() < capacity) {
				n.addTO(entry);
				if (n.getParent() != null) {
					System.out.println(KMAC1.ycDistance(n.getParent().getData(), entry.getData()) + " is the parent distance at leaf");
				}
			}
			else {
				this.split(n, entry);
			}
		}
		// keeps recursing into the leaf that is best for it
		else {
			ArrayList<TreeObject> objects = n.getObjects();
			double min = Double.MAX_VALUE;
			TreeObject mo = null; // this will be the object that it will be routed into
			System.out.println("running through all possiblities");
			for (TreeObject object: objects) {
				if (object.getChild() != null && object.getChild().getObjects().size() > 0) { // potential choice for mo
					double objDist = KMAC1.ycDistance(object.getData(), entry.getData());
					System.out.println(objDist);
					if (objDist < min) {
						min = objDist;
						mo = object;
					}
				}
			}
			System.out.println("ran through all possibilities");
			System.out.println(min + " was the parent distance at an ancestor");
			this.insertNode(mo.getChild(), entry);
		}
	}
	
	public void insertNode(MTreeNode n, TreeObject entry) {
		data.add(entry.getData().getIndex(), entry.getData());
		if (this.getRoot().isLeaf()) {
			this.insertNodeR(n, entry);
		}
		else {
			this.insertNodeNR(n, entry);
		}
	}
	
	public void insertNodeR(MTreeNode n, TreeObject entry) {
		entry.setContainer(n);
		if (n.getParent() != null) {
			System.out.println(KMAC1.ycDistance(n.getParent().getData(), entry.getData()) + " is the parent distance at leaf");
		}
		if (n.getObjects().size() < capacity) {
			n.addTO(entry);
		}
		else {
			this.split(n, entry);
		}
	}
	
	// leaves will have non-null parents guaranteed
	public void insertNodeNR(MTreeNode n, TreeObject entry) {
		data.add(entry.getData().getIndex(), entry.getData());
		double min = Double.MAX_VALUE;
		MTreeNode mo = null;
		for (MTreeNode leaf: n.getLeafNodes()) {
			double objDist = KMAC1.ycDistance(leaf.getParent().getData(), entry.getData());
			if (objDist < min) {
				min = objDist;
				mo = leaf;
			}
		}
		this.insertNodeR(mo, entry);
	}
	
	public int getSize() {
		return size;
	}
	
	public void setCurrent(int i) {
		currentEntry = i;
	}
	
	public int getCurrent() {
		return currentEntry;
	}
	
	class DistanceComparator implements Comparator<Pair<Integer, Double>> {
		@Override
		public int compare(Pair<Integer, Double> a, Pair<Integer, Double> b) {
			return a.getLast() < b.getLast() ? -1 : a.getLast().equals(b.getLast()) ? 0 : 1;
		}
	}
	
	// Makes a new object with children out of a list of TreeObjects. The returned object is one of the
	// members of the list (the "central" one), and the children are all the others.
	
	public TreeObject centralize(ArrayList<TreeObject> objects) {
		ArrayList<Double> distances = new ArrayList<Double>();
		for (int i = 0; i < objects.size(); i++) {
			distances.add(0.0);
		}
		for (int i = 0; i < objects.size(); i++) {
			for (int j = 0; j < i; j++) {
				double d = KMAC1.ycDistance(objects.get(i).getData(), objects.get(j).getData());
				distances.set(i, distances.get(i) + d);
				distances.set(j, distances.get(j) + d);
			}
		}
		int minIndex = 0;
		double minDist = distances.get(0);
		for (int i = 1; i < objects.size(); i++) {
			if (minDist > distances.get(i)) {
				minIndex = i;
				minDist = distances.get(i);
			}
		}
		TreeObject center = objects.get(minIndex);
		
		MTreeNode container = center.getContainer();
		objects.remove(minIndex);
		
		// center is now the central object from object, objects is the remaining ones
		
		if (center.getChild() == null || center.getChild().getObjects().size() == 0) { // need to make new child node
			MTreeNode newChild = new MTreeNode(true, center);
			for (TreeObject o: objects) {
				o.setContainer(newChild);
				if (o.getChild() != null) {
					if (o.getChild().getObjects().size() > 0) {
						newChild.setLeaf(false);
					}
				}
				newChild.getObjects().add(o);
			}
			center.setChild(newChild);
			container.clearTOs();
		}
		else if (center.getChild().getObjects().size() + objects.size() <= capacity) { // put objects in existing child node
			for (int i = 0; i < objects.size(); i++) {
				TreeObject o = objects.get(i);
				if (o.getChild() != null) {
					if (o.getChild().getObjects().size() > 0) {
						center.getChild().setLeaf(false);
					}
				}
				center.getChild().getObjects().add(o);
				o.setContainer(center.getChild());
			}
			container.clearTOs();
		}
		else { // centralize existing child node, put objects in existing child node which is 
			// guaranteed to not overflow because it has at most objects as the original node
			TreeObject returnedCenter = this.centralize(center.getChild().getObjects());
			objects.add(returnedCenter);
			MTreeNode newChild = new MTreeNode(true, center);
			center.setChild(newChild);
			newChild.setParent(center);
			for (TreeObject o: objects) {
				newChild.getObjects().add(o);
				o.setContainer(newChild);
				if (o.getChild() != null) {
					if (o.getChild().getObjects().size() > 0) {
						center.getChild().setLeaf(false);
					}
				}
			}
			container.clearTOs();
		}
		return center;
	}
	
	// Takes in the entries, "newEntries", of a node which needs to be split.
	// Returns promoted0 and promoted1, which contain as children the other entries in newEntries,
	// and resets all other entries in newEntries to have the correct parents.
	
	public Pair<TreeObject, TreeObject> splitHelper(ArrayList<TreeObject> newEntries) {
		// allocate new nodes partitionNode0, partitionNode1 with parents promoted0, promoted1
		
		Pair<TreeObject, TreeObject> promoted = MTree.optimalSplitPolicy(newEntries);
		
		TreeObject promoted0 = promoted.getFirst();
		ArrayList<TreeObject> promotedChildren0 = new ArrayList<TreeObject>();
		if (promoted0.getChild() != null) {
			promotedChildren0 = promoted0.getChild().getObjects();
		}
		TreeObject promoted1 = promoted.getLast();
		ArrayList<TreeObject> promotedChildren1 = new ArrayList<TreeObject>();
		if (promoted1.getChild() != null) {
			promotedChildren1 = promoted1.getChild().getObjects();
		}
		
		ArrayList<TreeObject> partition0 = new ArrayList<TreeObject>();
		ArrayList<TreeObject> partition1 = new ArrayList<TreeObject>();
		
		int numEntries = newEntries.size() - 2;
		
		LinkedList<Pair<Integer, Double>> distances0 = new LinkedList<Pair<Integer, Double>>();
		LinkedList<Pair<Integer, Double>> distances1 = new LinkedList<Pair<Integer, Double>>();
		
		// if n is a leaf or not
		for (int i = 0; i < newEntries.size(); i++) {
			TreeObject o = newEntries.get(i);
			double r = o.getR();
			if (o.getData().getIndex() != promoted0.getData().getIndex() && o.getData().getIndex() != promoted1.getData().getIndex()) {
				double d0 = KMAC1.ycDistance(promoted0.getData(), o.getData());
				distances0.add(new Pair<Integer, Double>(i, d0 + r));
				double d1 = KMAC1.ycDistance(promoted1.getData(), o.getData());
				distances1.add(new Pair<Integer, Double>(i, d1 + r));
			}
		}	
		
		Collections.sort(distances0, new DistanceComparator());
		Collections.sort(distances1, new DistanceComparator());
		
		for (int i = 0; i < numEntries - 2; i++) {
			if (i % 2 == 0) {
				int index = distances0.pop().getFirst();
				partition0.add(newEntries.get(index));
				for (int j = 0; j < distances1.size(); j++) {
					if (distances1.get(j).getFirst() == index) {
						distances1.remove(j);
						break;
					}
				}
			}
			else {
				int index = distances1.pop().getFirst();
				partition1.add(newEntries.get(index));
				for (int j = 0; j < distances0.size(); j++) {
					if (distances0.get(j).getFirst() == index) {
						distances0.remove(j);
						break;
					}
				}
			}
		}
		
		Pair<Integer, Double> biggest0 = distances0.pop();
		partition0.add(newEntries.get(biggest0.getFirst()));
		for (int j = 0; j < distances1.size(); j++) {
			if (distances1.get(j).getFirst() == biggest0.getFirst()) {
				distances1.remove(j);
				break;
			}
		}
		Pair<Integer, Double> biggest1 = distances1.pop();
		partition1.add(newEntries.get(biggest1.getFirst()));
		
		// make promoted0
		if (promotedChildren0.size() == 0) {
			MTreeNode partitionNode0 = new MTreeNode(true, promoted0);
			for (TreeObject o: partition0) {
				partitionNode0.addTO(o);
			}
			promoted0.setChild(partitionNode0);
			promoted0.getContainer().setLeaf(false);
		}
		else {
			if (promotedChildren0.size() + partition0.size() > capacity) {
				TreeObject center = this.centralize(promoted0.getChild().getObjects());
				// center is now all of promoted0's children into one node, promoted0's child currently 
				// is devoid of objects
				// what is center's current container? its old container, promoted0.getChild()
				partition0.add(center);
			}
			for (TreeObject o: partition0) {
				promoted0.getChild().addTO(o);
			}
			promoted0.getContainer().setLeaf(false);
		}
		
		// make promoted1
		if (promotedChildren1.size() == 0) {
			MTreeNode partitionNode1 = new MTreeNode(true, promoted1);
			for (TreeObject o: partition1) {
				partitionNode1.addTO(o);
			}
			promoted1.setChild(partitionNode1);
			promoted1.getContainer().setLeaf(false);
		}
		else {
			if (promotedChildren1.size() + partition1.size() > capacity) {
				TreeObject center = this.centralize(promoted1.getChild().getObjects());
				// center is now all of promoted1's children into one node, promoted1's child currently 
				// is devoid of objects
				// what is center's current container? its old container, promoted1.getChild()
				partition1.add(center);
			}
			for (TreeObject o: partition1) {
				promoted1.getChild().addTO(o);
			}
			promoted1.getContainer().setLeaf(false);
		}
		
		return new Pair<TreeObject, TreeObject>(promoted0, promoted1);
	}
	
	public void split(MTreeNode n, TreeObject entry) {
		//n.addTO(entry);
		ArrayList<TreeObject> newEntries = n.getObjects();
		newEntries.add(entry);
		entry.setContainer(n);
		
		TreeObject parentObject = null;
		MTreeNode parentContainer = null;
		TreeObject reAdd = null; // this is the TreeObject corresponding to parentObject if it exists, replaced by promoted0
		
		// if n is not the root
		
		if (n.getParent() != null) {
			parentObject = n.getParent();
			// need to make new object with parentObject's data to split
			reAdd = new TreeObject(n, null, parentObject.getData(), 0); 
			newEntries.add(reAdd);
			if (parentObject.getContainer() != null) { // removes parentObject from parentContainer
				parentContainer = parentObject.getContainer();
				for (int i = 0; i < parentContainer.getObjects().size(); i++) {
					if (parentContainer.getObjects().get(i).getData().getIndex() == (parentObject.getData().getIndex())) {
						parentContainer.getObjects().remove(i);
						break;
					}
				}
			}
		}

		Pair<TreeObject, TreeObject> promoted = this.splitHelper(newEntries);
		TreeObject promoted0 = promoted.getFirst(); // its container is still n
		TreeObject promoted1 = promoted.getLast(); // its container is still n
		
		if (n.getParent() == null) {
			MTreeNode newRoot = new MTreeNode(false, null);
			this.root = newRoot;
			newRoot.addTO(promoted0);
			newRoot.addTO(promoted1);
		}
		
		else {
			parentContainer.addTO(promoted0);
			if (parentContainer.getObjects().size() < capacity) {
				parentContainer.addTO(promoted1);
			}
			else {
				this.split(parentContainer, promoted1);
			}
		}
	}
	
	// fine, always returns 2 different indices with furthest possible distance
	public static Pair<TreeObject, TreeObject> optimalSplitPolicy(ArrayList<TreeObject> promote) {
		double furthestDist = 0;
		int k1 = promote.size();
		int k2 = promote.size();
		for (int i = 0; i < promote.size(); i++) {
			for (int j = 0; j < i; j++) {
				double ijDistance = KMAC1.ycDistance(promote.get(i).getData(), promote.get(j).getData());
				if (ijDistance > furthestDist) {
					k1 = i;
					k2 = j;
					furthestDist = ijDistance;
				}
			}
		}
		return new Pair<TreeObject, TreeObject>(promote.get(k1), promote.get(k2));
	}
	
	public Pair<Integer, ArrayList<Kmer>> rangeSearch(Kmer s, double r) {
		Pair<Integer, ArrayList<Kmer>> result = this.rangeSearch(s, r, root);
		/**
		try {
			String filename = "rangeSearchPruning.txt";
			FileWriter fw = new FileWriter(filename, true);
			fw.write(result.getFirst());
			fw.close();
		}
		catch (IOException ioe) {
			System.err.println("IOException: " + ioe.getMessage());
		}
		System.out.println("Range search has been completed at root.");
		System.out.println("----------------------------------------");
		System.out.println("Total number of entries pruned: " + result.getFirst());
		**/
		return result;
	}
	
	public Pair<Integer, ArrayList<Kmer>> rangeSearch(Kmer s, double r, MTreeNode n) {
		ArrayList<Kmer> inRange = new ArrayList<Kmer>();
		int pruned = 0;
		if (!n.isLeaf()) {
			for (TreeObject o: n.getObjects()) {
				if (KMAC1.ycDistance(o.getData(), s) <= r) {
					inRange.add(o.getData());
				}
				if (o.getChild() != null && o.getChild().getObjects().size() > 0) {
					if (KMAC1.ycDistance(o.getData(), s) <= r + o.getR())	{
						Pair<Integer, ArrayList<Kmer>> subRangeSearch = this.rangeSearch(s, r, o.getChild());
						inRange.addAll(subRangeSearch.getLast());
						pruned += subRangeSearch.getFirst();
					}
					else {
						//System.out.println("pruned successfully");
						pruned += o.getChild().getNumObjects();
					}
				}
			}
		}
		else {
			for (TreeObject o: n.getObjects()) {
				if (KMAC1.ycDistance(o.getData(), s) <= r) {
					inRange.add(o.getData());
				}
			}
			return new Pair<Integer, ArrayList<Kmer>>(pruned, inRange);
		}
		return new Pair<Integer, ArrayList<Kmer>>(pruned, inRange);
	}
	
	public ArrayList<Kmer> getData() {
		return data;
	}
	
	/**
	public static MTree constructTree(ArrayList<Kmer> kmers, int c, double pNR, int cu) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter("testDistribution.txt", "UTF-8");
		MTree tree = new MTree(c, pNR, cu, kmers.size());
		ArrayList<Kmer> kmerSet = new ArrayList<Kmer>();
		for (int i = 0; i < kmers.size(); i++) {
			System.out.println("beginning analysis of the " + i + " kmer");
			Kmer kmer = kmers.get(i);
			if (kmer instanceof GappedKmer) {
				System.out.println("gapped");
				for (Kmer k: ((GappedKmer) kmer).getBaseKmers()) {
					System.out.println(k);
				}
			}
			else {
				System.out.println("basic");
				System.out.println(kmer);
			}
			kmerSet.add(kmer);
		}
		for (int i = 0; i < kmers.size(); i++) {
			for (int j = 0; j < i; j++) {
				writer.println(KMAC1.ycDistance(kmerSet.get(i), kmerSet.get(j), pNR, cu));
			}
			System.out.println("making progress" + i);
		}
		System.out.println("reached");
		writer.close();
		return tree;
	}
	**/
	
	public static MTree constructTree(ArrayList<Kmer> kmers, int c) {
		MTree tree = new MTree(c, kmers.size());
		// for (int i = 0; i < kmers.size(); i++) {
		for (int i = 0; i < kmers.size(); i++) {
			tree.setCurrent(i);
			System.out.println("start"+i);
			Kmer data = kmers.get(i);
			data.setIndex(i);
			tree.insertNode(tree.getRoot(), tree.new TreeObject(null, null, data, 0));
			System.out.println(tree.root.getNumObjects());
		}
		System.out.println(kmers.size());
		System.out.println(tree.root.getNumObjects());
		for (TreeObject leaf: tree.root.getLeaves()) {
			leaf.setParentsR();
		}
		tree.root.recursiveSetNumObjects();
		// tree.root.recursiveRI();
		// tree.root.recursivePrint();
		return tree;
	}
	
	public void setSize(int s) {
		this.size = s;
	}
	
	public void setRoot(MTreeNode n) {
		this.root = n;
	}
	
	// pre: s1.length() == s2.length()
	public static int test1(String s1, String s2) {
		int mm = 0;
		for (int i = 0; i < s1.length(); i++) {
			if (s1.charAt(i) != s2.charAt(i)) {
				mm++;
			}
		}
		return mm;
	}
	
	public static int test2(String s1, String s2) {
		if (s1.length() > s2.length()) {
			return MTree.test2(s2, s1);
		}
		else {
			int mm = s2.length();
			for (int i = 1 - s1.length(); i < s2.length(); i++) {
				int start2 = Math.max(0, i);
				int end2 = Math.min(i + s1.length() - 1, s2.length() - 1);
				int start1 = -1;
				int end1 = -1;
				if (start2 == 0) {
					start1 = s1.length() - 1 - (end2 - start2);
					end1 = s1.length() - 1;
				}
				else if (end2 == s2.length() - 1) {
					start1 = 0;
					end1 = end2 - start2;
				}
				else {
					start1 = 0;
					end1 = s1.length() - 1;
				}
				String ss1 = s1.substring(start1, end1 + 1);
				String ss2 = s2.substring(start2, end2 + 1);
				mm = Math.min(mm, s2.length() - (end2 - start2 + 1) 
						+ s1.length() - (end1 - start1 + 1) + MTree.test1(ss1, ss2));
			}
			return mm;
		}
	}
	
	public static int testDist(String s1, String s2) {
		s1 = s1.toUpperCase();
		s2 = s2.toUpperCase();
		String s1r = SequenceUtils.reverseComplement(s1);
		// return Math.min(Math.min(MTree.test2(s1, s2, cutoff), MTree.test2(s1r, s2, cutoff)), cutoff);
		return Math.min(MTree.test2(s1, s2), MTree.test2(s1r, s2));
	}
	
	public static int testDistWC(String s1, String s2, int cutoff) {
		// this one has cutoff
		s1 = s1.toUpperCase();
		s2 = s2.toUpperCase();
		String s1r = SequenceUtils.reverseComplement(s1);
		return Math.min(Math.min(MTree.test2(s1, s2), MTree.test2(s1r, s2)), cutoff);
	}
	
	public static void main(String[] args) {
		System.out.println(MTree.test1("hey", "bey"));
		System.out.println(MTree.test1("hey", "bcy"));
		System.out.println(MTree.test2("substringABCD", "ABCsubs"));
	}
}
