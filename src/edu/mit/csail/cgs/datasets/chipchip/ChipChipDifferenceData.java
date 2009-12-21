/*
 * Created on Nov 28, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.utils.NotFoundException;

public class ChipChipDifferenceData implements ChipChipData {
    
    private ChipChipData d1, d2;
    
    public ChipChipDifferenceData(ChipChipData c1, ChipChipData c2) { 
        d1 = c1; 
        d2 = c2;
    }

    public int getExptID(int i, int j) {
        return -1;
    }

    public double getIP(int i, int j) {
        return d1.getIP(i, j) - d2.getIP(i, j);
    }

    public double getMax() {
    	double max = (-1*Double.MAX_VALUE);
    	for (int i = 0; i < getCount(); i++) {
            for (int j = 0; j < getReplicates(i); j++) {
            	double curr = getRatio(i,j);
            	if(curr>max){
            		max =curr; 
            	}
            }
    	}return(max);
    }
    public double getMin() {
    	double min = Double.MAX_VALUE;
    	for (int i = 0; i < getCount(); i++) {
            for (int j = 0; j < getReplicates(i); j++) {
            	double curr = getRatio(i,j);
            	if(curr<min){
            		min =curr; 
            	}
            }
    	}return(min);
    }

    public double getMax(String chrom, int start, int stop) throws NotFoundException {
    	double max = (-1*Double.MAX_VALUE);
    	for (int i = 0; i < getCount(); i++) {
            for (int j = 0; j < getReplicates(i); j++) {
            	double curr = getRatio(i,j);
            	if(curr>max){
            		max =curr; 
            	}
            }
    	}return(max);
    }  
    public double getMin(String chrom, int start, int stop) throws NotFoundException {
    	double min = Double.MAX_VALUE;
    	for (int i = 0; i < getCount(); i++) {
            for (int j = 0; j < getReplicates(i); j++) {
            	double curr = getRatio(i,j);
            	if(curr<min){
            		min =curr; 
            	}
            }
    	}return(min);
    }
    public double getRatio(int i, int j) {
        return d1.getRatio(i, j) - d2.getRatio(i, j);
    }

    public char getStrand(int i, int j) {
        return d1.getStrand(i, j);
    }

    public double getVar(int i, int j) {
        return d1.getVar(i, j) - d2.getVar(i, j);
    }

    public double getWCE(int i, int j) {
        return d1.getWCE(i, j) - d2.getWCE(i, j);
    }

    public void window(String chrom, int start, int stop, double minratio) throws NotFoundException {
        d1.window(chrom, start, stop);
        d2.window(chrom, start, stop);
    }

    public int baseToIndex(int b) {
        return d1.baseToIndex(b);
    }

    public int getCount() {
        return d1.getCount();
    }

    public String getName() {
        return String.format("Diff(%s,%s)", d1.getName(), d2.getName());
    }

    public int getPos(int i) {
        return d1.getPos(i);
    }

    public int getReplicates(int i) {
        return d1.getReplicates(i);
    }

    public double getValue(int i, int j) {
        return d1.getValue(i, j) - d2.getValue(i, j);
    }

    public String getVersion() {
        return String.format("%s:%s", d1.getVersion(), d2.getVersion());
    }

    public String getWindowChrom() {
        return d1.getWindowChrom();
    }

    public int getWindowStart() {
        return d1.getWindowStart();
    }

    public int getWindowStop() {
        return d1.getWindowStop();
    }

    public void window(String chrom, int start, int stop) throws NotFoundException {
        d1.window(chrom, start, stop);
        d2.window(chrom, start, stop);
    }

}
