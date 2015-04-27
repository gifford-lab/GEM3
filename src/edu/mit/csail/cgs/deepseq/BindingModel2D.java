package edu.mit.csail.cgs.deepseq;

import java.util.List;
import java.util.ArrayList;

import edu.mit.csail.cgs.utils.Pair;

public class BindingModel2D {

    protected int min, max;     // the start and end position
    protected int summit;       // the position of highest prob point
    
    protected ArrayList<Pair<Integer, List<Pair<Integer, Double>>>> empiricalDistribution;
    
    public BindingModel2D(ArrayList<Pair<Integer, List<Pair<Integer, Double>>>> bindingDist){
        min=0; max=0;
        empiricalDistribution=bindingDist;
        //loadData(bindingDist);
        //makeProbabilities();
    }
}
