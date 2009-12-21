package edu.mit.csail.cgs.ewok.verbs.expression;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;

import edu.mit.csail.cgs.datasets.expression.ExpressionLoader;
import edu.mit.csail.cgs.datasets.expression.LocatedProbe;
import edu.mit.csail.cgs.datasets.expression.ProbePlatform;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.utils.Closeable;
import edu.mit.csail.cgs.utils.iterators.EmptyIterator;

public class ExpressionRegionExpander 
	implements Expander<Region,LocatedProbe>, Closeable {
	
	private boolean shouldCloseLoader;
	private ExpressionLoader loader;
	private ProbePlatform platform;

	public ExpressionRegionExpander(String pp) throws SQLException {
		shouldCloseLoader = true;
		loader = new ExpressionLoader();
		platform = loader.loadPlatform(pp);
	}
	
	public ExpressionRegionExpander(ExpressionLoader el, String pp) throws SQLException { 
		shouldCloseLoader = false;
		loader = el;
		platform = loader.loadPlatform(pp);
	}
	
	public ExpressionLoader getLoader() { return loader; }

	public Iterator<LocatedProbe> execute(Region a) {
		try {
			Collection<LocatedProbe> list = loader.loadProbesInRegion(a, platform);
			return list.iterator();
			
		} catch (SQLException e) {
			e.printStackTrace();
			return new EmptyIterator<LocatedProbe>();
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
