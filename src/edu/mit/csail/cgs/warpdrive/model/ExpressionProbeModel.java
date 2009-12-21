/**
 * 
 */
package edu.mit.csail.cgs.warpdrive.model;

import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;

/**
 * @author tdanford
 */
public class ExpressionProbeModel extends RegionExpanderModel<LocatedExprMeasurement> {
    public ExpressionProbeModel(Expander<Region,LocatedExprMeasurement> exp) { 
        super(exp);
    }
}
