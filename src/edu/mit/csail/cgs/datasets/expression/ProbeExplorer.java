package edu.mit.csail.cgs.datasets.expression;

import java.sql.*;
import java.util.*;

/**
 * 
 * @author tdanford
 * @see Experiment
 * @see ExpressionLoader
 * @see ExprMeasurement
 * @see Probe
 */
public class ProbeExplorer {
	
	public static void main(String[] args) { 
		String name = args.length > 0 ? args[0] : "Hoxa1";
		try {
			ProbeExplorer pe = new ProbeExplorer(new ExpressionLoader());
			pe.exploreProbe(name);
			/*
			pe.addExperiment("PPG_relative_tp1_bg0");
			ExpressionGenerator eeg = pe.createGenerator();
			Set<String> names = new HashSet<String>(); names.add("Hoxa1"); 
			eeg.getValues(names);
			*/
			pe.close();
			
		
			
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

	private ExpressionLoader loader;
	private Vector<Experiment> expts;

	public ProbeExplorer(ExpressionLoader l) { 
		loader = l;
		expts = null;
	}
	
	/**
	 * adds an <tt>Experiment</tt> with name <tt>name</tt> to <tt>expts</tt> field which is a <tt>Vector</tt> of <tt>Experiments</tt>.
	 * @param name name of the experiment. E.g. one name could be "<tt>PPG_relative_tp1_bg0</tt>"
	 * @throws SQLException
	 */
	public void addExperiment(String name) throws SQLException { 
		Experiment e = loader.loadExperiment(name);
		addExperiment(e);
	}

	/**
	 * adds an <tt>Experiment e</tt> to <tt>expts</tt> field which is a <tt>Vector</tt> of <tt>Experiments</tt>.
	 * @param e the experiment as an <tt>Experiment</tt> object
	 */
	public void addExperiment(Experiment e) { 
		if(expts == null) { expts = new Vector<Experiment>(); }
		if(!expts.contains(e)) { expts.add(e); }
	}
	
	/**
	 * closes the <tt>ProbeExplorer</tt>. Necessary to end connection with the database.
	 */
	public void close() { 
		loader.close();
	}
	
	/**
	 * Gets all the probes that have the name <tt>n</tt>. For each of these probes, 
	 * it finds the experiments that each probe participates in and outputs the name 
	 * of each experiment and the (expression) value of the probe (in this experiment) 
	 * on the console.
	 * @param n name of the probe
	 * @throws SQLException
	 */
	public void exploreProbe(String n) throws SQLException { 
		Collection<Probe> probes = loader.loadProbesByName(n);
		System.out.println("Found " + probes.size() + " probes.");
		for(Probe p : probes) { 
			System.out.println(p.getName() + " (#" + p.getDBID() + ")");
			Collection<Experiment> pexpts = loader.loadAllExperiments(p);
			for(Experiment expt : pexpts) { 
				if(expts == null || expts.contains(expt)) { 
					ExprMeasurement m = loader.loadMeasurement(expt, p);
					double value = m.getValue();
					System.out.println("\t" + expt.getName() + " --> " + value);
				}
			}
		}
	}
	
	/**
	 * creates a generator which implements ExpressionGenerator
	 * @return
	 */
	public ExpressionGenerator createGenerator() { 
		return new ExperimentExpressionGenerator();
	}

	public static interface ExpressionGenerator { 
		public Map<String,Map<String,Double>> getValues(Set<String> names);
	}

	private class ExperimentExpressionGenerator implements ExpressionGenerator { 

		public ExperimentExpressionGenerator() {}
		
		/**
		 * 
		 * @param names Names of the probes as a <tt>Set</tt>
		 * @return a map. Each entry of the map contains the name of the probe as a key and all the experiments that
		 * this probe participates in and its corresponding values as a value.
		 */
		public Map<String,Map<String,Double>> getValues(Set<String> names) { 
			TreeMap<String,Map<String,Double>> map = new TreeMap<String,Map<String,Double>>();
			if(expts == null) { return map; }

			for(String n : names) { 
				try { 
   					Collection<Probe> probes = loader.loadProbesByName(n);
					if(probes.size() > 0) { 
						map.put(n, new TreeMap<String,Double>());
						Probe p = probes.iterator().next();

						for(Experiment expt : expts) { 
							ExprMeasurement m = loader.loadMeasurement(expt, p);
							if(m != null) { 
								double value = m.getValue();
								map.get(n).put(expt.getName(),value);
							} else { 
								map.get(n).put(expt.getName(), 0.0);
							}
						}

					}
				} catch(SQLException se) { 
					se.printStackTrace(System.err); 
					if(map.containsKey(n)) { map.remove(n); }
				}
			}
			return map;
		}
	}
}
