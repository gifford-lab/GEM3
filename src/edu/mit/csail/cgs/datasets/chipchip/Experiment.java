package edu.mit.csail.cgs.datasets.chipchip;

public class Experiment {
    protected String name, version, replicate;
    protected int dbid;
    protected int fragdist, factorone, factortwo, cellsone, cellstwo, conditionone, conditiontwo;
    protected int species;
    protected boolean active;

    public String getName() {return name;}
    public String getVersion() {return version;}
    public String getReplicate() {return replicate;}
    public int getDBID() {return dbid;}
    public int getSpecies() {return species;}
    public int getFragDist() {return fragdist;}
    public int getFactorOne() {return factorone;}
    public int getCellsOne() {return cellsone;}
    public int getConditionOne() {return conditionone;}
    public int getFactorTwo() {return factortwo;}
    public int getCellsTwo() {return cellstwo;}
    public int getConditionTwo() {return conditiontwo;}
    public boolean getActive() {return active;}

}
