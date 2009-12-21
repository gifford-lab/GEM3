package edu.mit.csail.cgs.utils;

import java.util.*;

/* from a set of values for positive and negative examples, computes 
   ROC curve areas, returns plots, and the sensitivity and specificity values
*/

public class ROC {

    private ArrayList<SS> vals;

    /* creates a new ROC object with a set of positive and negative examples.
       The double values are whatever value the threshold can vary over */
    public ROC (double[] positive, double[] negative, Comparator c) {
        double min, max, step;
        min = positive[0]; max = positive[0];
        for (int i = 0; i < positive.length; i++) {
            if (positive[i] > max) {max = positive[i];}
            if (positive[i] < min) {min = positive[i];}
        }
        for (int i = 0; i < negative.length; i++) {
            if (negative[i] > max) {max = negative[i];}
            if (negative[i] < min) {min = negative[i];}
        }
        step = (max - min) / 100;
        vals = new ArrayList<SS>();
        int i = 0;
        for (double threshold = min - step; threshold <= max + step; threshold += step) {
            int tp = 0, tn = 0, fp = 0, fn = 0; 
            for (int j = 0; j < positive.length; j++) {
                if (c.compare(positive[j],threshold) == 1) {
                    tp++;
                } else {
                    fn++;
                }
            }
            for (int j = 0; j < negative.length; j++) {
                if (c.compare(negative[j],threshold) == 1){
                    fp++;
                } else {
                    tn++;
                }
            }
            SS n = new SS();
            vals.add(n);
            n.sens = tp / (tp+fn+.000001);
            n.spec = tn / (tn+fp+.000001);
            i++;
        }
    }    

    /* 2D ROC.  second dimension must have size 2.
       Takes two comparators.  The comparator will be called on 
       the value and the threshold and should return 1 iff the
       value meets the threshold */
    public ROC(double[][] positive, double[][] negative, Comparator c1, Comparator c2) {
        double minone, mintwo, maxone, maxtwo;
        minone = positive[0][0];
        mintwo = positive[0][1];
        maxone = positive[0][0];
        maxtwo = positive[0][1];
        for (int i = 0; i < positive.length; i++) {
            if (positive[i][0] > maxone) {maxone = positive[i][0];}
            if (positive[i][0] < minone) {minone = positive[i][0];}
            if (positive[i][1] > maxtwo) {maxtwo = positive[i][1];}
            if (positive[i][1] < mintwo) {mintwo = positive[i][1];}
        }
        for (int i = 0; i < negative.length; i++) {
            if (negative[i][0] > maxone) {maxone = negative[i][0];}
            if (negative[i][0] < minone) {minone = negative[i][0];}
            if (negative[i][1] > maxtwo) {maxtwo = negative[i][1];}
            if (negative[i][1] < mintwo) {mintwo = negative[i][1];}
        }
        double step1 = (maxone - minone)/20;
        double step2 = (maxtwo - mintwo)/20;
        if (step2 == 0) {
            step2 = .01;
            step1 = (maxone - minone) / 400;
        }
        vals = new ArrayList<SS>();
        System.err.println("Minone=" +minone +" maxone=" + maxone + " step1="+ step1);
        System.err.println("Mintwo=" +mintwo +" maxtwo=" + maxtwo + " step2="+ step2);
        for (double t1 = minone; t1 <= maxone; t1 += step1) {
            for (double t2 = mintwo; t2 <= maxtwo; t2 += step2) {
                int tp = 0, tn = 0, fp = 0, fn = 0; 
                for (int j = 0; j < positive.length; j++) {
                    if (c1.compare(positive[j][0],t1) == 1 &&
                        c2.compare(positive[j][1],t2) == 1) {
                        tp++;
                    } else {
                        fn++;
                    }
                }
                for (int j = 0; j < negative.length; j++) {
                    if (c1.compare(negative[j][0],t1) == 1 &&
                        c2.compare(negative[j][1],t2) == 1) {
                        fp++;
                    } else {
                        tn++;
                    }
                }
                SS n = new SS();
                vals.add(n);
                n.sens = tp / (tp + fn + .000001);
                n.spec = tn / (tn + fp + .000001);
                //                System.err.println("At t1=" + t1 + " and t2=" + t2 + "  sens=" + n.sens + " and spec="+ n.spec);
            }
        }
    }

    public double getROCArea () {
        double area = 0;
        SS vals[] = new SS[this.vals.size()];
        this.vals.toArray(vals);
        Arrays.sort(vals);
        for (int i = 0 ; i < vals.length - 1; i++) {
            area += .5 * (vals[i].spec + vals[i+1].spec) * (vals[i+1].sens - vals[i].sens);
        }
        if (vals[0].spec > vals[vals.length-1].spec) {
            area += vals[0].sens * .5 * ( 1 + vals[0].spec);
            area += (1 - vals[vals.length-1].sens) * .5 * vals[vals.length-1].spec;
        } else {
            area += vals[vals.length-1].sens * .5 * ( 1 + vals[vals.length-1].spec);
            area += (1 - vals[0].sens) * .5 * vals[0].spec;
        }
        if (area < .5) {
            return 1 - area;
        } else {
            return area;
        }
    }

}

class SS implements Comparable{
    public double sens, spec;
    public int compareTo(Object o) {
        double val = ((SS)o).sens;
        if (sens < val) {
            return -1;
        } else if (sens > val) {
            return 1;
        } else {
            return 0;
        }
    }
}

