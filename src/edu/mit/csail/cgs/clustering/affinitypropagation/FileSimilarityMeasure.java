package edu.mit.csail.cgs.clustering.affinitypropagation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;

import edu.mit.csail.cgs.clustering.Clusterable;
import edu.mit.csail.cgs.clustering.SimpleClusterable;

/**
 * @author reeder
 *
 */
public class FileSimilarityMeasure<X extends Clusterable> extends SimilarityMeasure<X> {

	Vector<Vector<Double>> simValues = new Vector<Vector<Double>>();
	Vector<Clusterable> objects = new Vector<Clusterable>();
	HashMap<Pair, Double> valuemap;
	double prefvalue;
	
	public FileSimilarityMeasure(String simfile, String divider, double prefvalue) {
		this.prefvalue = prefvalue;
		String line = "init";
		String[] splitline = new String[1];
		valuemap = new HashMap<Pair, Double>();
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(simfile));
			while (!((line = dataIn.readLine())==null)) {
				splitline = line.split(divider);
				SimpleClusterable object1 = new SimpleClusterable(splitline[0]);
				SimpleClusterable object2 = new SimpleClusterable(splitline[1]);
				if (!objects.contains(object1)) {
					objects.add(object1);
				}
				if (!objects.contains(object2)) {
					objects.add(object2);
				}
				valuemap.put(new Pair(object1, object2), Double.valueOf(splitline[2]));
			}
			System.err.println("Measures: "+valuemap.size());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	/*
	public FileSimilarityMeasure(String simfile, String preffile, String divider, int notused, int notused2) {
		String line = "init";
		String[] splitLine = new String[1];
		int objcnt = 0;
		int idx0, idx1;
		HashMap<String, Integer> idxmap = new HashMap<String, Integer>();
		try {
			BufferedReader dataIn = new BufferedReader(new FileReader(simfile));
			while (!((line = dataIn.readLine())==null)) {
				splitLine = line.split(divider);
				if (!idxmap.containsKey(splitLine[0])) {
					idxmap.put(splitLine[0], objcnt++);
					names.add(splitLine[0]);
				}
				if (!idxmap.containsKey(splitLine[1])) {
					idxmap.put(splitLine[1], objcnt++);
					names.add(splitLine[1]);
				}
				idx0 = idxmap.get(splitLine[0]);
				idx1 = idxmap.get(splitLine[1]);
				if (simValues.size()<=idx0) {
					simValues.setSize(idx0+1);
				}
				if (simValues.get(idx0)==null) {
					simValues.set(idx0, new Vector<Double>());
				}
				if (simValues.get(idx0).size()<=idx1) {
					simValues.get(idx0).setSize(idx1+1);
				}
				simValues.get(idx0).set(idx1, Double.valueOf(splitLine[2]));
			}
			simValues.get(simValues.size()-1).setSize(simValues.size());
			BufferedReader prefin = new BufferedReader(new FileReader(preffile));
			idx0 = 0;
			while (!((line = prefin.readLine())==null)) {
				simValues.get(idx0).set(idx0++, Double.valueOf(line));
			}
		} catch (Exception e) {
			System.out.println(line);
			System.out.println(splitLine[0]);
			e.printStackTrace();
		}
	}
	*/
	
	@Override
	public void addNoise() {
		noiseAdded = true;

	}

	@Override
	public double get(int idx0, int idx1) {
		System.out.println("get shouldn't be called");
		if (simValues.get(idx0).get(idx1)==null) {
			return NEGINF;
		} else {
			return simValues.get(idx0).get(idx1);
		}
	}

	@Override
	public int size() {
		return objects.size();
	}
	
	public String getName(int idx) {
		return objects.get(idx).name();
	}

	/*public double evaluate(Clusterable e1, Clusterable e2) {
		if (valuemap.containsKey(e1.name()+e2.name()));
		return 0;
	}*/

	public double evaluate(X e1, X e2) {
		if (e1.name().equals(e2.name())) {
			return prefvalue;
		} else if (valuemap.containsKey(e1.name()+"SEP"+e2.name())) {
			return valuemap.get(e1.name()+"SEP"+e2.name());
		} else if (valuemap.containsKey(e1.name()+"SEP"+e1.name())) {
			return valuemap.get(e1.name()+"SEP"+e2.name());
		} else {
			return NEGINF;
		}
	}
	
	public double evaluate(Pair p) {
		if (p.symmetric()) {
			return prefvalue;
		} else if (valuemap.containsKey(p)) {
			return valuemap.get(p);
		} else {
			return NEGINF;
		}
	}
	
	public boolean exists(X e1, X e2) {
		return (e1.name().equals(e2.name())) || (valuemap.containsKey(e1.name()+"SEP"+e2.name())) || 
			(valuemap.containsKey(e1.name()+"SEP"+e1.name()));
	}
	
	public boolean exists(Pair p) {
		return p.symmetric() || valuemap.containsKey(p);
	}
	
	public Vector<Clusterable> objects() {
		return objects;
	}
	
}
