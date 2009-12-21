/*
 * Created on Mar 20, 2006
 */
package edu.mit.csail.cgs.datasets.chipchip;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.io.*;

import edu.mit.csail.cgs.datasets.species.Genome;

/**
 * @author tdanford
 */
public class MSPProbe extends Probe {
    
    private static NumberFormat nf;
    
    static { 
        nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(2);
    }

    /**
     * @param genome
     * @param chrom
     * @param loc
     */
    public MSPProbe(Genome genome, String chrom, int loc, String key, ChipChipMSP msp, int index) {
        super(genome, chrom, loc);
        addMSPValue(key, msp, index);
    }
    
    public void addMSPValue(String k, ChipChipMSP msp, int index) { 
        double[] array = new double[5];
        array[0] = msp.getRatio(index);
        array[1] = msp.getX(index);
        array[2] = msp.getPval(index);
        array[3] = msp.getPval3(index);
        array[4] = msp.getMedianOfRatios(index);
        super.addValue(k, array);
    }
    
    public double getRatio(String k) { return getValue(k)[0]; }
    public double getX(String k) { return getValue(k)[1]; }
    public double getPval(String k) { return getValue(k)[2]; }
    public double getPval3(String k) { return getValue(k)[3]; }
    public double getMoR(String k) { return getValue(k)[4]; }
    
    public double getMaxMoR() { 
        double val = 0.0;
        for(String k : getKeys()) { 
            if(getMoR(k) > val) { val = getMoR(k); }
        }
        return val;
    }
    
    public void debugPrint(PrintStream ps) {
        NumberFormat nf = DecimalFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        
        ps.println("Probe (" + getLocationString() + ")");
        for(String key : getKeys()) { 
            ps.print("\t" + key + " --> ");
            ps.print(" RoM: " + nf.format(getRatio(key)));
            ps.print(" MoR: " + nf.format(getRatio(key)));
            ps.print(" X: " + nf.format(getRatio(key)));
            ps.print(" Pval: " + nf.format(getRatio(key)));
            ps.print(" Pval3: " + nf.format(getRatio(key)));
            ps.println();
        }
    }
    
    public String toString() { 
        StringBuilder sb = new StringBuilder();
        sb.append(super.toString());
        for(String key : getKeys()) { 
            sb.append(" (" + key + ": ");
            sb.append(nf.format(getMoR(key)));
            sb.append(" ");
            sb.append(nf.format(getPval3(key)));
            sb.append(",");
            sb.append(nf.format(getPval(key)));
            sb.append(")");
        }
        return sb.toString();
    }
}
