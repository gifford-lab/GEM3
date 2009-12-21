package edu.mit.csail.cgs.ewok.verbs.expression;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;

import edu.mit.csail.cgs.datasets.expression.Experiment;
import edu.mit.csail.cgs.datasets.expression.ExprMeasurement;
import edu.mit.csail.cgs.datasets.expression.ExpressionLoader;
import edu.mit.csail.cgs.datasets.expression.LocatedExprMeasurement;
import edu.mit.csail.cgs.datasets.expression.LocatedProbe;
import edu.mit.csail.cgs.datasets.expression.ProbePlatform;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class LocatedExprMeasurementExpander 
	implements Expander<Region,LocatedExprMeasurement>, Closeable {
	
	private ExpressionLoader loader;
	private boolean shouldCloseLoader;
    private Experiment expt;
    private ProbePlatform plat;
	
	public LocatedExprMeasurementExpander(ExpressionLoader l, 
			String exptName, String platName) throws SQLException { 
		shouldCloseLoader = false;
		loader = l;
        expt = loader.loadExperiment(exptName);
        plat = loader.loadPlatform(platName);
	}

	public Iterator<LocatedExprMeasurement> execute(Region a) {        
        try {
            Collection<LocatedExprMeasurement> exprs = loader.loadMeasurementsInRegion(a, plat, expt);
            return exprs.iterator();        
        } catch (SQLException e) {
            e.printStackTrace();
            return new EmptyIterator<LocatedExprMeasurement>();
        }
	}

	public void close() {
		if(shouldCloseLoader && !loader.isClosed()) { 
			loader.close();
		}
		loader.close();
	}

	public boolean isClosed() {
		return loader == null || loader.isClosed();
	}
}
