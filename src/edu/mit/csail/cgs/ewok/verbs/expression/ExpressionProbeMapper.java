package edu.mit.csail.cgs.ewok.verbs.expression;

import java.sql.SQLException;

import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.datasets.expression.ExprMeasurement;
import edu.mit.csail.cgs.datasets.expression.ExpressionLoader;
import edu.mit.csail.cgs.datasets.expression.Probe;
import edu.mit.csail.cgs.ewok.verbs.Filter;
import edu.mit.csail.cgs.utils.Closeable;

public class ExpressionProbeMapper 
	implements Filter<Probe,ExprMeasurement>, Closeable {
	
	private ExpressionLoader loader;
	private boolean shouldCloseLoader;
	private Experiment expt;
	
	public ExpressionProbeMapper(ExpressionLoader el, String exptName) throws SQLException { 
		loader = el;
		shouldCloseLoader = false;
		expt = loader.loadExperiment(exptName);
	}
	
	public ExpressionProbeMapper(String exptName) throws SQLException { 
		loader = new ExpressionLoader();
		shouldCloseLoader = true;
		expt = loader.loadExperiment(exptName);
	}
	
	public ExpressionLoader getLoader() { return loader; }
	public Experiment getExperiment() { return expt; }

	public ExprMeasurement execute(Probe a) {
		try {
			return loader.loadMeasurement(expt, a);
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}
	}

	public void close() {
		if(shouldCloseLoader && !loader.isClosed()) { 
			loader.close();
		}
		loader = null;
	}

	public boolean isClosed() {
		return loader == null || loader.isClosed();
	}
}
