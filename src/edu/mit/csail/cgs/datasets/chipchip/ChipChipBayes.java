package edu.mit.csail.cgs.datasets.chipchip;

import edu.mit.csail.cgs.utils.NotFoundException;

public interface ChipChipBayes extends GenericExperiment {    
    public void window(String chrom, int start, int stop, double minpost, double minsize) throws NotFoundException;
    public double getPosterior(int i);
    public double getStrength(int i);
    public double getPosteriorStd(int i);
    public double getStrengthStd(int i);    
    public double getMaxStrength();
    public double getMaxStrength(String chrom, int start, int stop) throws NotFoundException;
    public double getMaxPosterior();
    public double getMaxPosterior(String chrom, int start, int stop) throws NotFoundException;
}
