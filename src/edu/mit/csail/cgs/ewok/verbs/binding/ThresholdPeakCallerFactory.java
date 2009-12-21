/*
 * Created on Sep 29, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.binding;

import edu.mit.csail.cgs.datasets.chipchip.Probe;
import edu.mit.csail.cgs.datasets.locators.ChipChipLocator;
import edu.mit.csail.cgs.datasets.locators.ExptLocator;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.PeakCallerFactory;
import edu.mit.csail.cgs.ewok.verbs.probers.ChipChipImmediateProbeGenerator;

/**
 * @author tdanford
 */
public class ThresholdPeakCallerFactory implements PeakCallerFactory {
    
    public ThresholdPeakCallerFactory() {
    }

    /* (non-Javadoc)
     * @see edu.mit.csail.cgs.ewok.PeakCallerFactory#createCaller(edu.mit.csail.cgs.Bio.Genome, edu.mit.csail.cgs.datasets.locators.ExptLocator)
     */
    public PeakCaller createCaller(Genome g, ExptLocator loc, Object args) {
        ChipChipLocator aloc = (ChipChipLocator)loc;
        ChipChipImmediateProbeGenerator gen = new ChipChipImmediateProbeGenerator(g, aloc);
        RegionProber<Probe> prober = new RegionProber.Wrapper<Probe>(gen);
        
        String ename = aloc.name + "," + aloc.version;
        double thresh = args != null ? (Double)args : 2.0;
        ThresholdPeakFinder<Probe> finder = new ThresholdPeakFinder<Probe>(ename, thresh);
        
        PeakCaller caller = new PeakCaller.FromFinder<Probe>(prober, finder);
        return caller;
    }

}
