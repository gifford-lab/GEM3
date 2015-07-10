package edu.mit.csail.cgs.deepseq.discovery.kmer.mtree;

import java.io.BufferedWriter;
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


public class MTree {
	
	private MTreeNode root;
	private int capacity;
	private double posNegRatio;
	private int cutoff;
	private int size;
	private ArrayList<Kmer> data;
	
	class MTreeNode {
		
		private boolean isLeaf;
		private ArrayList<TreeObject> nro;
		private TreeObject parent;
		
		public MTreeNode(boolean l, TreeObject p) {
			isLeaf = l;
			parent = p;
			nro = new ArrayList<TreeObject>();
		}
		
		public MTreeNode(boolean l, TreeObject p, ArrayList<TreeObject> objects) {
			isLeaf = l;
			parent = p;
			nro = objects;
		}
		
		public void addTO(TreeObject to) {
			if (to.getContainer() != null) {
				MTreeNode container = to.getContainer();
				if (container.getParent() != null) {
					System.out.println(container.getParent().getData().getLast() + " is the old container");
				}
				int containerSize = container.getObjects().size();
				for (int i = 0; i < containerSize; i++) {
					if (container.getObjects().get(i).getData().getLast().equals(to.getData().getLast())) {
						container.getObjects().remove(i);
						System.out.println("successful");
						break;
					}
				}
			}
			nro.add(to);
			to.setContainer(this);
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
		
		public ArrayList<Pair<Kmer, Integer>> getData() {
			ArrayList<Pair<Kmer, Integer>> data = new ArrayList<Pair<Kmer, Integer>>();
			for (TreeObject o: this.nro) {
				data.add(o.getData());
			}
			return data;
		}
		
		public int getNumObjects() {
			int total = this.nro.size();
			if (!this.isLeaf()) {
				for (TreeObject o: this.nro) {
					if (o.getChild() != null) {
						total += o.getChild().getNumObjects();
					}
				}
			}
			return total;
		}
		
		public void recursiveRI() {
			if (this.getParent() == null) {
				System.out.println("BEGINNING RECURSIVE RI ---------------");
			}
			for (TreeObject o: this.getObjects()) {
				if (o.getChild() != null) {
					/**
					if (this.getParent() != null) {
						if (o.getData().getLast().equals(this.getParent().getData().getLast())) {
							break;
						}
					}
					**/
					o.getChild().recursiveRI();
				}
				if (this.getParent() == null) {
					System.out.println(o.getData().getLast() + " parent: none");
				}
				else {
					System.out.println(o.getData().getLast() + " parent: " + this.getParent().getData().getLast());
				}
				o.RI();
			}
		}
		
		public void splitRecursiveRI() {
			if (this.getParent() == null) {
				System.out.println("BEGINNING SPLIT RECURSIVE RI ---------------");
			}
			for (TreeObject o: this.getObjects()) {
				if (o.getChild() != null) {
					o.getChild().recursiveRI();
				}
				if (this.getParent() == null) {
					System.out.println(o.getData().getLast() + " parent: none");
				}
				else {
					System.out.println(o.getData().getLast() + " parent: " + this.getParent().getData().getLast());
				}
				o.RI();
			}
		}
	}
	
	class TreeObject {
		
		private MTreeNode container;
		private double radius; // covering radius
		private MTreeNode child; // covering tree
		private Pair<Kmer, Integer> data;
		private double pnRatio;
		private int cut;
		private double parentDistance;
		
		public TreeObject(MTreeNode co, MTreeNode c, Pair<Kmer, Integer> d, double r, double pnr, int cu) {
			container = co;
			child = c;
			data = d;
			pnRatio = pnr;
			cut = cu;
			if (co == null) {
				parentDistance = Integer.MAX_VALUE;
			}
			else if (co.getParent() != null) {
				parentDistance = KMAC1.ycDistance(co.getParent().getData().getFirst(), d.getFirst(), pnRatio, cut);
			}
			else {
				parentDistance = Integer.MAX_VALUE;
			}
			radius = r;
		}
		
		public TreeObject(MTreeNode co, Pair<Kmer, Integer> d) {
			container = co;
			child = null;
			data = d;
			pnRatio = 1;
			cut = Integer.MAX_VALUE;
			if (co.getParent() != null) {
				parentDistance = KMAC1.ycDistance(co.getParent().getData().getFirst(), d.getFirst(), pnRatio, cut);
			}
			else {
				parentDistance = Integer.MAX_VALUE;
			}
			radius = 0;
		}
		
		public void RI() {
			if (this.child != null) {
				if (this.child.getObjects().size() == 7) {
					System.out.println("FIRST INSTANCE OF BAD");
				}
			}
		}
		
		public Pair<Kmer, Integer> getData() {
			return data;
		}
		
		public void resetParentDistance() {
			if (this.container.getParent() != null) {
				parentDistance = KMAC1.ycDistance(this.container.getParent().getData().getFirst(), this.data.getFirst(), this.pnRatio, this.cut);
			}
		}
		
		public void setData(Pair<Kmer, Integer> d) {
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
		
		public double getParentDistance() {
			return parentDistance;
		}
		
		public MTreeNode getContainer() {
			return container;
		}
		
		public void setContainer(MTreeNode newContainer) {
			container = newContainer;
		}
		
		public double getPNR() {
			return this.pnRatio;
		}
		
		public int getCutoff() {
			return this.cut;
		}
	}
	
	public MTree(int k, double pNR, int c, int s) {
		capacity = k;
		posNegRatio = pNR;
		cutoff = c;
		root = new MTreeNode(true, null);
		size = s;
		data = new ArrayList<Kmer>();
		for (int i = 0; i < size; i++) {
			data.add(new Kmer());
		}
	}
	
	public MTreeNode getRoot() {
		return root;
	}
	
	public void insertNode(MTreeNode n, TreeObject entry) {
		System.out.println("inserting id number " + entry.getData().getLast() + " into node child of " + (n.getParent() != null ? n.getParent().getData().getLast() : "null"));
		this.root.recursiveRI();
		System.out.println("insertion");
		data.add(entry.getData().getLast(), entry.getData().getFirst());
		if (n.isLeaf()) {
			entry.setContainer(n);
			System.out.println("leafinsertion");
			if (n.getObjects().size() < capacity) {
				System.out.println("fine");
				n.addTO(entry);
				entry.resetParentDistance();
				entry.RI();
			}
			else {
				if (n.getParent() != null) {
					System.out.println("entering leaf of " + n.getParent().getData().getLast());
				}
				System.out.println("leafsplit");				
				Pair<TreeObject, TreeObject> splitResult = this.split(n, entry);
				entry.RI();
			}
		}
		else {
			System.out.println("childfinder");
			ArrayList<TreeObject> objects = n.getObjects();
			ArrayList<TreeObject> containers = new ArrayList<TreeObject>();
			double min = Double.MAX_VALUE;
			TreeObject mo = null;
			for (TreeObject object: objects) {
				double objDist = KMAC1.ycDistance(object.getData().getFirst(), entry.getData().getFirst(), this.posNegRatio, this.cutoff);
				if (objDist < min) {
					min = objDist;
					mo = object;
				}
				if (objDist <= object.getR()) {
					containers.add(object);
				}
			}
			if (containers.isEmpty()) {
				mo.setR(KMAC1.ycDistance(mo.getData().getFirst(), entry.getData().getFirst(), this.posNegRatio, this.cutoff));
			}
			System.out.println("TYPE 1 inserting in " + mo.getData().getLast());
			this.insertNode(mo.getChild(), entry);
			entry.RI();
		}
	}
	
	public int getSize() {
		return size;
	}
	
	class DistanceComparator implements Comparator<Pair<Integer, Double>> {
		@Override
		public int compare(Pair<Integer, Double> a, Pair<Integer, Double> b) {
			return a.getLast() < b.getLast() ? -1 : a.getLast().equals(b.getLast()) ? 0 : 1;
		}
	}
	
	public ArrayList<TreeObject> splitHelper(ArrayList<TreeObject> entries) {
		ArrayList<TreeObject> splitResult = new ArrayList<TreeObject>();
		
		return splitResult;
	}
	
	public Pair<TreeObject, TreeObject> multiSplit(MTreeNode n, ArrayList<TreeObject> entries) {
		this.root.splitRecursiveRI();
		System.out.println("multisplit");
		
		//n.addTO(entry);
		ArrayList<TreeObject> newEntries = n.getObjects();
		for (TreeObject o: entries) {
			n.addTO(o);
		}
		
		TreeObject parentObject = null;
		MTreeNode parentContainer = null;
		TreeObject reAdd = null; // this is the TreeObject corresponding to parentObject if it exists, replaced by promoted0
		
		// if n is not the root
		
		if (n.getParent() != null) {
			parentObject = n.getParent();
			if (parentObject.getContainer() != null) {
				parentContainer = parentObject.getContainer();
			}
			reAdd = new TreeObject(null, null, parentObject.getData(), 0, this.posNegRatio, this.cutoff);
		}
		
		// allocate new nodes partitionNode0, partitionNode1 with parents promoted0, promoted1
		
		ArrayList<TreeObject> promoted = MTree.optimalSplitPolicy(newEntries, this.posNegRatio, this.cutoff);
		TreeObject promoted0 = promoted.get(0);
		ArrayList<TreeObject> promotedChildren0 = new ArrayList<TreeObject>();
		if (promoted0.getChild() != null) {
			promotedChildren0 = promoted0.getChild().getObjects();
		}
		TreeObject promoted1 = promoted.get(1);
		ArrayList<TreeObject> promotedChildren1 = new ArrayList<TreeObject>();
		if (promoted1.getChild() != null) {
			promotedChildren1 = promoted1.getChild().getObjects();
		}
		
		// System.out.println(promotedChildren0.size() + promotedChildren1.size() + "GGGGGG");
		
		double maxR0 = 0;
		double maxR1 = 0;
		
		ArrayList<TreeObject> partition0 = new ArrayList<TreeObject>();
		ArrayList<TreeObject> partition1 = new ArrayList<TreeObject>();
		
		int numEntries = newEntries.size() - 2;
		
		LinkedList<Pair<Integer, Double>> distances0 = new LinkedList<Pair<Integer, Double>>();
		LinkedList<Pair<Integer, Double>> distances1 = new LinkedList<Pair<Integer, Double>>();
		
		// if n is a leaf or not
		
		if (n.isLeaf()) {
			for (int i = 0; i < newEntries.size(); i++) {
				TreeObject o = newEntries.get(i);
				if (o.getData().getLast() != promoted0.getData().getLast() && o.getData().getLast() != promoted1.getData().getLast()) {
					double d0 = KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances0.add(new Pair<Integer, Double>(i, d0));
					double d1 = KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances1.add(new Pair<Integer, Double>(i, d1));
				}
			}			
		}
		else {
			for (int i = 0; i < newEntries.size(); i++) {
				TreeObject o = newEntries.get(i);
				double r = o.getR();
				if (o.getData().getLast() != promoted0.getData().getLast() && o.getData().getLast() != promoted1.getData().getLast()) {
					double d0 = KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances0.add(new Pair<Integer, Double>(i, d0 + r));
					double d1 = KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances1.add(new Pair<Integer, Double>(i, d1 + r));
				}
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
					}
				}
			}
			else {
				int index = distances1.pop().getFirst();
				partition1.add(newEntries.get(index));
				for (int j = 0; j < distances0.size(); j++) {
					if (distances0.get(j).getFirst() == index) {
						distances0.remove(j);
					}
				}
			}
		}
		
		Pair<Integer, Double> biggest0 = distances0.pop();
		partition0.add(newEntries.get(biggest0.getFirst()));
		for (TreeObject pentry: partition0) {
			System.out.println(pentry.getData().getLast());
		}
		maxR0 = biggest0.getLast();
		System.out.println("not broken " + promoted0.getData().getLast());
		for (int j = 0; j < distances1.size(); j++) {
			if (distances1.get(j).getFirst() == biggest0.getFirst()) {
				distances1.remove(j);
			}
		}
		Pair<Integer, Double> biggest1 = distances1.pop();
		partition1.add(newEntries.get(biggest1.getFirst()));
		for (TreeObject pentry: partition1) {
			System.out.println(pentry.getData().getLast());
		}
		maxR1 = biggest1.getLast();
		System.out.println("not broken pair " + promoted1.getData().getLast());
		
		// at this point all elements have been added to partition0 and partition1
		/**
		for (TreeObject o: partition0) {
			o.setContainer(null);
		}
		
		for (TreeObject o: partition1) {
			o.setContainer(null);
		}
		*/
		
		// n.getObjects().clear();
		
		if (promotedChildren0.size() == 0) {
			System.out.println("XX0");
			MTreeNode partitionNode0 = new MTreeNode(n.isLeaf(), promoted0);
			for (TreeObject o: partition0) {
				partitionNode0.addTO(o);
				o.resetParentDistance();
			}
			promoted0.setR(maxR0);
			promoted0.setChild(partitionNode0);
			promoted0.RI();
		}
		else {
			System.out.println("YY0");
			int alreadyHave = promoted0.getChild().getObjects().size();
			if (alreadyHave + partition0.size() <= capacity) {
				for (int i = 0; i < partition0.size(); i++) {
					promoted0.getChild().addTO(partition0.get(i));
					System.out.println("normalinsertion0");
				}
				maxR0 = 0;
				for (TreeObject o: promoted0.getChild().getObjects()) {
					maxR0 = Math.max(KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff) + o.getR(), maxR0);
				}
				promoted0.setR(maxR0);
				promoted0.RI();
			}
			else {
				this.multiSplit(promoted0.getChild(), partition0);
				System.out.println("weirdinsertion0");
			}
		}
		
		if (promotedChildren1.size() == 0) {
			System.out.println("XX1");
			MTreeNode partitionNode1 = new MTreeNode(n.isLeaf(), promoted1);
			for (TreeObject o: partition1) {
				partitionNode1.addTO(o);
				o.resetParentDistance();
			}
			promoted1.setR(maxR1);
			promoted1.setChild(partitionNode1);
			promoted1.RI();
		}
		else {
			System.out.println("YY1");
			int alreadyHave = promoted1.getChild().getObjects().size();
			System.out.println("WE ALREADY HAVE " + alreadyHave);
			if (alreadyHave + partition1.size() <= capacity) {
				for (int i = 0; i < partition1.size(); i++) {
					promoted1.getChild().addTO(partition1.get(i));
					System.out.println("normalinsertion1");
				}
				maxR1 = 0;
				for (TreeObject o: promoted1.getChild().getObjects()) {
					System.out.println("which part?1");
					maxR1 = Math.max(KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff) + o.getR(), maxR0);
				}
				promoted1.setR(maxR1);
				promoted1.RI();
			}
			else {
				Pair<TreeObject, TreeObject> promotedPair = this.multiSplit(promoted1.getChild(), partition1);
				System.out.println("weirdinsertion1");
			}
		}
		
		if (n.getParent() == null) {
			MTreeNode newRoot = new MTreeNode(false, null);
			this.root = newRoot;
			newRoot.addTO(promoted0);
			promoted0.resetParentDistance();
			newRoot.addTO(promoted1);
			promoted1.resetParentDistance();
			
			System.out.println("ROOT WAS RESET.");
		}
		
		else {
			//System.out.println("drop");
			if (parentContainer == null) {
				System.out.println("break0");
			}
			for (int i = 0; i < parentContainer.getObjects().size(); i++) {
				if (parentContainer.getObjects().get(i).getData().getLast().equals(parentObject.getData().getLast())) {
					parentContainer.getObjects().remove(i);
				}
			}
			parentContainer.addTO(promoted0);
			if (promoted1.getContainer() != null) {
				MTreeNode container = promoted1.getContainer();
				int containerSize = container.getObjects().size();
				for (int i = 0; i < containerSize; i++) {
					if (container.getObjects().get(i).getData().getLast().equals(promoted1.getData().getLast())) {
						container.getObjects().remove(i);
						promoted1.setContainer(null);
						break;
					}
				}
			}
			
			promoted0.resetParentDistance();
			if (parentContainer.getObjects().size() < capacity) {
				parentContainer.addTO(promoted1);
				promoted1.resetParentDistance();
			}
			else {
				this.split(parentContainer, promoted1);
			}
			System.out.println("readd");
			System.out.println("TYPE 4");
			this.insertNode(this.root, reAdd);
		}
		
		return new Pair<TreeObject, TreeObject>(promoted0, promoted1);
	}
	
	public Pair<TreeObject, TreeObject> split(MTreeNode n, TreeObject entry) {
		this.root.splitRecursiveRI();
		System.out.println("split");
		
		//n.addTO(entry);
		ArrayList<TreeObject> newEntries = n.getObjects();
		newEntries.add(entry);
		
		TreeObject parentObject = null;
		MTreeNode parentContainer = null;
		TreeObject reAdd = null; // this is the TreeObject corresponding to parentObject if it exists, replaced by promoted0
		
		// if n is not the root
		
		if (n.getParent() != null) {
			parentObject = n.getParent();
			if (parentObject.getContainer() != null) {
				parentContainer = parentObject.getContainer();
			}
			reAdd = new TreeObject(null, null, parentObject.getData(), 0, this.posNegRatio, this.cutoff);
		}
		
		// allocate new nodes partitionNode0, partitionNode1 with parents promoted0, promoted1
		
		ArrayList<TreeObject> promoted = MTree.optimalSplitPolicy(newEntries, this.posNegRatio, this.cutoff);
		TreeObject promoted0 = promoted.get(0);
		ArrayList<TreeObject> promotedChildren0 = new ArrayList<TreeObject>();
		if (promoted0.getChild() != null) {
			promotedChildren0 = promoted0.getChild().getObjects();
		}
		TreeObject promoted1 = promoted.get(1);
		ArrayList<TreeObject> promotedChildren1 = new ArrayList<TreeObject>();
		if (promoted1.getChild() != null) {
			promotedChildren1 = promoted1.getChild().getObjects();
		}
		
		// System.out.println(promotedChildren0.size() + promotedChildren1.size() + "GGGGGG");
		
		double maxR0 = 0;
		double maxR1 = 0;
		
		ArrayList<TreeObject> partition0 = new ArrayList<TreeObject>();
		ArrayList<TreeObject> partition1 = new ArrayList<TreeObject>();
		
		int numEntries = newEntries.size() - 2;
		
		LinkedList<Pair<Integer, Double>> distances0 = new LinkedList<Pair<Integer, Double>>();
		LinkedList<Pair<Integer, Double>> distances1 = new LinkedList<Pair<Integer, Double>>();
		
		// if n is a leaf or not
		
		if (n.isLeaf()) {
			for (int i = 0; i < newEntries.size(); i++) {
				TreeObject o = newEntries.get(i);
				if (o.getData().getLast() != promoted0.getData().getLast() && o.getData().getLast() != promoted1.getData().getLast()) {
					double d0 = KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances0.add(new Pair<Integer, Double>(i, d0));
					double d1 = KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances1.add(new Pair<Integer, Double>(i, d1));
				}
			}			
		}
		else {
			for (int i = 0; i < newEntries.size(); i++) {
				TreeObject o = newEntries.get(i);
				double r = o.getR();
				if (o.getData().getLast() != promoted0.getData().getLast() && o.getData().getLast() != promoted1.getData().getLast()) {
					double d0 = KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances0.add(new Pair<Integer, Double>(i, d0 + r));
					double d1 = KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff);
					distances1.add(new Pair<Integer, Double>(i, d1 + r));
				}
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
					}
				}
			}
			else {
				int index = distances1.pop().getFirst();
				partition1.add(newEntries.get(index));
				for (int j = 0; j < distances0.size(); j++) {
					if (distances0.get(j).getFirst() == index) {
						distances0.remove(j);
					}
				}
			}
		}
		
		Pair<Integer, Double> biggest0 = distances0.pop();
		partition0.add(newEntries.get(biggest0.getFirst()));
		for (TreeObject pentry: partition0) {
			System.out.println(pentry.getData().getLast());
		}
		maxR0 = biggest0.getLast();
		System.out.println("not broken " + promoted0.getData().getLast());
		for (int j = 0; j < distances1.size(); j++) {
			if (distances1.get(j).getFirst() == biggest0.getFirst()) {
				distances1.remove(j);
			}
		}
		Pair<Integer, Double> biggest1 = distances1.pop();
		partition1.add(newEntries.get(biggest1.getFirst()));
		for (TreeObject pentry: partition1) {
			System.out.println(pentry.getData().getLast());
		}
		maxR1 = biggest1.getLast();
		System.out.println("not broken pair " + promoted1.getData().getLast());
		
		// at this point all elements have been added to partition0 and partition1
		/**
		for (TreeObject o: partition0) {
			o.setContainer(null);
		}
		
		for (TreeObject o: partition1) {
			o.setContainer(null);
		}
		*/
		
		// n.getObjects().clear();
		
		if (promotedChildren0.size() == 0) {
			System.out.println("XX0");
			MTreeNode partitionNode0 = new MTreeNode(n.isLeaf(), promoted0);
			for (TreeObject o: partition0) {
				partitionNode0.addTO(o);
				o.resetParentDistance();
			}
			promoted0.setR(maxR0);
			promoted0.setChild(partitionNode0);
			promoted0.RI();
		}
		else {
			System.out.println("YY0");
			int alreadyHave = promoted0.getChild().getObjects().size();
			System.out.println("WE ALREADY HAVE " + alreadyHave);
			for (int i = 0; i < Math.min(capacity - alreadyHave, partition0.size()); i++) {
				promoted0.getChild().addTO(partition0.get(i));
				System.out.println("normalinsertion0");
			}
			for (int i = capacity - alreadyHave; i < partition0.size() - capacity + alreadyHave; i++) {
				System.out.println("TYPE 2");
				// trying out new stuff
				TreeObject addition = partition0.get(i);
				for (int j = 0; j < newEntries.size(); j++) {
					if (newEntries.get(j).getData().getLast().equals(addition.getData().getLast())) {
						newEntries.remove(j);
					}
				}
				addition.setContainer(null);
				this.insertNode(promoted0.getChild(), addition);
				// this.insertNode(promoted0.getChild(), partition0.get(i));
				System.out.println("weirdinsertion0");
			}
			/**
			for (TreeObject o: partition0) {
				this.insertNode(promoted0.getChild(), o);
			}
			**/
			maxR0 = 0;
			for (TreeObject o: promoted0.getChild().getObjects()) {
				System.out.println("which part?0");
				maxR0 = Math.max(KMAC1.ycDistance(promoted0.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff) + o.getR(), maxR0);
			}
			promoted0.setR(maxR0);
			promoted0.RI();
		}
		
		System.out.println(partition0.size()+"AAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		System.out.println(partition1.size()+"BBBBBBBBBBBBBBBBBBBBBBBBBBBB");
		
		if (promotedChildren1.size() == 0) {
			System.out.println("XX1");
			MTreeNode partitionNode1 = new MTreeNode(n.isLeaf(), promoted1);
			for (TreeObject o: partition1) {
				partitionNode1.addTO(o);
				o.resetParentDistance();
			}
			promoted1.setR(maxR1);
			promoted1.setChild(partitionNode1);
			promoted1.RI();
		}
		else {
			System.out.println("YY1");
			int alreadyHave = promoted1.getChild().getObjects().size();
			for (int i = 0; i < Math.min(capacity - alreadyHave, partition1.size()); i++) {
				promoted1.getChild().addTO(partition1.get(i));
				System.out.println("normalinsertion1");
			}
			for (int i = capacity - alreadyHave; i < partition1.size() - capacity + alreadyHave; i++) {
				System.out.println("TYPE 3");
				// trying out new stuff
				TreeObject addition = partition1.get(i);
				for (int j = 0; j < newEntries.size(); j++) {
					if (newEntries.get(j).getData().getLast().equals(addition.getData().getLast())) {
						newEntries.remove(j);
					}
				}
				addition.setContainer(null);
				this.insertNode(promoted1.getChild(), addition);
				// this.insertNode(promoted1.getChild(), partition1.get(i));
				System.out.println("weirdinsertion1");
			}
			/**
			for (TreeObject o: partition1) {
				this.insertNode(promoted1.getChild(), o);
			}
			**/
			maxR1 = 0;
			for (TreeObject o: promoted1.getChild().getObjects()) {
				maxR1 = Math.max(KMAC1.ycDistance(promoted1.getData().getFirst(), o.getData().getFirst(), this.posNegRatio, this.cutoff) + o.getR(), maxR1);
			}
			promoted1.setR(maxR1);
			promoted1.RI();
		}
		
		if (n.getParent() == null) {
			MTreeNode newRoot = new MTreeNode(false, null);
			this.root = newRoot;
			newRoot.addTO(promoted0);
			promoted0.resetParentDistance();
			newRoot.addTO(promoted1);
			promoted1.resetParentDistance();
			
			System.out.println("ROOT WAS RESET.");
		}
		
		else {
			//System.out.println("drop");
			if (parentContainer == null) {
				System.out.println("break0");
			}
			for (int i = 0; i < parentContainer.getObjects().size(); i++) {
				if (parentContainer.getObjects().get(i).getData().getLast().equals(parentObject.getData().getLast())) {
					parentContainer.getObjects().remove(i);
				}
			}
			parentContainer.addTO(promoted0);
			if (promoted1.getContainer() != null) {
				MTreeNode container = promoted1.getContainer();
				int containerSize = container.getObjects().size();
				for (int i = 0; i < containerSize; i++) {
					if (container.getObjects().get(i).getData().getLast().equals(promoted1.getData().getLast())) {
						container.getObjects().remove(i);
						promoted1.setContainer(null);
						break;
					}
				}
			}
			
			promoted0.resetParentDistance();
			if (parentContainer.getObjects().size() < capacity) {
				parentContainer.addTO(promoted1);
				promoted1.resetParentDistance();
			}
			else {
				System.out.println("also");
				if (parentContainer.getObjects().size() != 6) {
					System.out.println("why oh why " + parentContainer.getObjects().size());
				}
				this.split(parentContainer, promoted1);
				System.out.println("also this");
			}
			System.out.println("readd");
			System.out.println("TYPE 4");
			this.insertNode(this.root, reAdd);
		}
		
		return new Pair<TreeObject, TreeObject>(promoted0, promoted1);
	}
	
	// fine, always returns 2 different indices with furtheste possible distance
	public static ArrayList<TreeObject> optimalSplitPolicy(ArrayList<TreeObject> promote, double pNR, int c) {
		ArrayList<TreeObject> promoted = new ArrayList<TreeObject>();
		double furthestDist = 0;
		int k1 = promote.size();
		int k2 = promote.size();
		for (int i = 0; i < promote.size(); i++) {
			for (int j = 0; j < i; j++) {
				double ijDistance = KMAC1.ycDistance(promote.get(i).getData().getFirst(), promote.get(j).getData().getFirst(), pNR, c);
				if (ijDistance > furthestDist) {
					k1 = i;
					k2 = j;
					furthestDist = ijDistance;
				}
			}
		}
		promoted.add(promote.get(k1));
		promoted.add(promote.get(k2));
		return promoted;
	}
	
	public ArrayList<Pair<Kmer, Integer>> rangeSearch(Kmer s, double r) {
		return this.rangeSearch(s, r, root);
	}
	
	public ArrayList<Pair<Kmer, Integer>> rangeSearch(Kmer s, double r, MTreeNode n) {
		ArrayList<Pair<Kmer, Integer>> inRange = new ArrayList<Pair<Kmer, Integer>>();
		if (!n.isLeaf()) {
			for (TreeObject o: n.getObjects()) {
				double objDist = KMAC1.ycDistance(o.getData().getFirst(), s, this.posNegRatio, this.cutoff);
				if (objDist <= r) {
					inRange.add(o.getData());
				}
				if (Math.abs(o.getParentDistance() - KMAC1.ycDistance(n.getParent().getData().getFirst(), s, this.posNegRatio, this.cutoff)) <= r + o.getR()) {
					double subDist = KMAC1.ycDistance(o.getData().getFirst(), s, this.posNegRatio, this.cutoff);
					if (subDist <= r + o.getR()) {
						inRange.addAll(this.rangeSearch(s, r, o.getChild()));
					}
				}
			}
		}
		else {
			for (TreeObject o: n.getObjects()) {
				if (Math.abs(o.getParentDistance() - KMAC1.ycDistance(n.getParent().getData().getFirst(), s, this.posNegRatio, this.cutoff)) <= r) {
					double subDist = KMAC1.ycDistance(o.getData().getFirst(), s, this.posNegRatio, this.cutoff);
					if (subDist <= r) {
						inRange.add(o.getData());
					}
				}
			}
			return inRange;
		}
		return inRange;
	}
	
	public ArrayList<Kmer> getData() {
		return data;
	}
	
	public static MTree constructTree(ArrayList<Kmer> kmers, int c, double pNR, int cu) {
		MTree tree = new MTree(c, pNR, cu, kmers.size());
		for (int i = 0; i < kmers.size(); i++) {
			System.out.println("start"+i);
			Pair<Kmer, Integer> data = new Pair<Kmer, Integer>(kmers.get(i), i); 
			tree.insertNode(tree.getRoot(), tree.new TreeObject(null, null, data, 0, pNR, cu));
			System.out.println(tree.root.getNumObjects());
		}
		System.out.println(kmers.size());
		System.out.println(tree.root.getNumObjects());
		return tree;
	}
	
	public void setSize(int s) {
		this.size = s;
	}
	
	public void setRoot(MTreeNode n) {
		this.root = n;
	}
	
	public static void main(String[] args) {
		ArrayList<Kmer> kmers = new ArrayList<Kmer>();
	}
}
