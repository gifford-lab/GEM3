package edu.mit.csail.cgs.warpdrive.model;

import java.util.*;
import edu.mit.csail.cgs.utils.*;
import edu.mit.csail.cgs.datasets.chipseq.*;
import edu.mit.csail.cgs.datasets.general.Region;

public class ChipSeqAnalysisModel extends WarpModel implements RegionModel, Runnable {

    private ChipSeqAnalysis analysis;
    private ChipSeqAnalysisProperties props;

    private Collection<ChipSeqAnalysisResult> results;
    private Region region;
    private boolean newinput;

    public ChipSeqAnalysisModel(ChipSeqAnalysis a) {
        analysis = a;
        props = new ChipSeqAnalysisProperties();
    }
    public void clearValues() {
        results = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public boolean isReady() {return !newinput;}
    public Collection<ChipSeqAnalysisResult> getResults() { return results;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    results = analysis.getResults(region);
                } catch (Exception e) {
                    e.printStackTrace();
                    results = new ArrayList<ChipSeqAnalysisResult>();
                }
            }
        }
    }

}