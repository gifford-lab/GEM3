package edu.mit.csail.cgs.datasets.chipchip;

public class FragDist {
    protected String name, version, description;
    protected int dbid;
    protected double[] values;

    public String getName() {return name;}
    public String getVersion() {return version;}
    public String getDescription() {return description;}
    public int getDBID() {return dbid;}
    public int getLength() {return values.length;}
    public double getValue(int i) {return values[i];}
}
