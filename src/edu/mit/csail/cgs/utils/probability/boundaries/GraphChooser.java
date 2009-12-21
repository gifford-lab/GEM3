/*
 * Created on Feb 10, 2006
 */
package edu.mit.csail.cgs.utils.probability.boundaries;

import java.util.*;

/**
 * @author tdanford
 */
public class GraphChooser implements ConstrainedChooser {
    
    private CountingGraph graph;

    public GraphChooser() {
        graph = new CountingGraph(6, 3);
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.psrg.tdanford.boundary.pvalues.ConstrainedChooser#logConstrainedChoose(int, int, int, int, boolean)
     */
    public double logConstrainedChoose(int N, int E, int p0, int s0, boolean optimal) {
        graph = new CountingGraph(N, E);
        graph.arrange(N, E, p0, s0, !optimal);
        return graph.getLogCount();
    }

}
