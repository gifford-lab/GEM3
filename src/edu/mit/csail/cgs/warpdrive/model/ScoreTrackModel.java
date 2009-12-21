/**
 * 
 */
package edu.mit.csail.cgs.warpdrive.model;

import edu.mit.csail.cgs.datasets.chippet.ChipPetDatum;
import edu.mit.csail.cgs.datasets.expression.*;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.ewok.verbs.Expander;
import edu.mit.csail.cgs.ewok.verbs.Mapper;
import edu.mit.csail.cgs.utils.Interval;

/**
 * @author tdanford
 */
public class ScoreTrackModel extends RegionExpanderModel<ChipPetDatum> {
    public ScoreTrackModel(Expander<Region,ChipPetDatum> exp) { 
        super(exp);
    }
}
