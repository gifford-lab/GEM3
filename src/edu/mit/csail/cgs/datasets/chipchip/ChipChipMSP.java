package edu.mit.csail.cgs.datasets.chipchip;

import java.util.*;
import edu.mit.csail.cgs.utils.parsing.*;
import edu.mit.csail.cgs.utils.*;

public interface ChipChipMSP extends GenericExperiment {
    public float getRatio(int i);
    public float getX(int i);
    public float getPval(int i);
    public float getPval3(int i);
    public float getMedianOfRatios(int i);


}
