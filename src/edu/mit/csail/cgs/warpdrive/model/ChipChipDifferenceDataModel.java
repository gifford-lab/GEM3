package edu.mit.csail.cgs.warpdrive.model;

import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import edu.mit.csail.cgs.utils.NotFoundException;
import edu.mit.csail.cgs.utils.Pair;
import edu.mit.csail.cgs.warpdrive.paintable.WarpPaintable;
import edu.mit.csail.cgs.datasets.chipchip.ChipChipData;
import edu.mit.csail.cgs.datasets.chipchip.GenericExperiment;
import edu.mit.csail.cgs.datasets.general.Region;

public class ChipChipDifferenceDataModel extends WarpModel implements Runnable, RegionModel {

    private ChipChipData data1, data2;
    private boolean newregion;
    private Region region;
    
    private Vector<Pair<Integer,Double>> diffs;
    private double max;

    public ChipChipDifferenceDataModel(ChipChipData d1, ChipChipData d2) {
        super();
        data1 = d1;
        data2 = d2;
        max = 0.0;
        newregion = false;
        
        diffs = new Vector<Pair<Integer,Double>>();
    }
    
    public Vector<Pair<Integer,Double>> getDiffs() { 
        return new Vector<Pair<Integer,Double>>(diffs);
    }

    public double getMax() {
        return max;
    }
    
    public synchronized void run() {
        while (keepRunning()) {
            try {
                if (!newregion) {
                    wait();
                }
            } catch (InterruptedException ex) {
                
            }
            if (newregion) {
                try {
                    max = 0.0;
                    
                    data1.window(region.getChrom(),region.getStart(),region.getEnd());
                    data2.window(region.getChrom(),region.getStart(),region.getEnd());
                    
                    diffs.clear();
                    
                    if(data1.getCount() == data2.getCount()) { 

                        Vector<Integer> pos = new Vector<Integer>();
                        Vector<Double> vals1 = new Vector<Double>();
                        Vector<Double> vals2 = new Vector<Double>();
                        
                        double s1 = 0.0, s2 = 0.0;
                        double logMaxMult = 0.0;
                        int c = data1.getCount();

                        for(int i = 0; i < c; i++) {
                            int p1 = data1.getPos(i);
                            int p2 = data2.getPos(i);

                            if(p1 == p2) {
                                Pair<Double,Double> av1 = getAveragePair(data1, i);
                                Pair<Double,Double> av2 = getAveragePair(data2, i);
                                double v1 = av1.getLast(), v2 = av2.getLast();
                                
                                double ratio = Math.log(v2 > 0.0 ? v1 / v2 : 1.0);
                                if(Math.abs(ratio) > Math.abs(logMaxMult)) { logMaxMult = ratio; }
                                
                                s1 += v1; s2 += v2;
                                vals1.add(v1); 
                                vals2.add(v2);
                                pos.add(p1);
                            }
                        }
                        
                        s1 /= (c > 0 ? (double)c : 1.0);
                        s2 /= (c > 0 ? (double)c : 1.0);

                        //double v2Mult = s2 > 0.0 ? s1 / s2 : 1.0;
                        //double v2Mult = Math.exp(logMaxMult);
                        double v2Mult = 1.0;
                        
                        for(int i = 0; i < pos.size(); i++) { 

                            int ploc = pos.get(i);
                            double v1 = vals1.get(i); 
                            double v2 = vals2.get(i) * v2Mult;
                            double diff = v1 - v2;
                            max = Math.max(max, Math.abs(diff));
                               
                            Pair<Integer,Double> p = new Pair<Integer,Double>(ploc, diff);
                            diffs.add(p);
                        }
                    }
                    
                    newregion = false;
                    notifyListeners();                    
                } catch (NotFoundException ex) {
                    // don't do anything here.  If we got an invalid chromosome, nothing
                    // will happen and the user will notice (hopefully)
                    ex.printStackTrace();
                }
            }
        }
    }
    
    private synchronized double getAverageValue(GenericExperiment expt, int index) { 
        double sum = 0.0;
        int repls = expt.getReplicates(index);
        for(int i = 0; i < repls; i++) { 
            sum += expt.getValue(index, i);
        }
        sum /= (double)repls;
        return sum;
    }

    private synchronized double getAverageIP(ChipChipData expt, int index) { 
        double sum = 0.0;
        int repls = expt.getReplicates(index);
        for(int i = 0; i < repls; i++) { 
            sum += expt.getIP(index, i);
        }
        sum /= (double)repls;
        return sum;
    }

    private synchronized Pair<Double,Double> getAveragePair(ChipChipData expt, int index) { 
        double ratiosum = 0.0;
        double ipsum = 0.0;
        int repls = expt.getReplicates(index);
        for(int i = 0; i < repls; i++) { 
            ipsum += expt.getIP(index, i);
            ratiosum += expt.getRatio(index, i);
        }
        ipsum /= (double)repls;
        ratiosum /= (double)repls;
        return new Pair<Double,Double>(ipsum, ratiosum);
    }
    
    public synchronized void setRegion(Region r) {
        if (newregion == false) {
            newregion = true;
            region = r;
        }
    }
    
    public Region getRegion() { return region; }
    public boolean isReady() { return !newregion; }
    
}
