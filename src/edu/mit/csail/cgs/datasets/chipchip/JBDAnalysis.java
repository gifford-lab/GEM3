package edu.mit.csail.cgs.datasets.chipchip;

import java.util.Map;
import java.util.Set;
import edu.mit.csail.cgs.datasets.species.Organism;

public class JBDAnalysis extends NameVersion {
    protected int dbid;
    protected Set<Experiment> inputs;
    protected Map<String,String> params;
    protected Organism species;

    public JBDAnalysis(Organism org, String name, String version) {
        super(name,version);
        dbid = -1;
        species = org;
    }
    public JBDAnalysis(Organism org, String name, String version, int dbid) {
        super(name,version);
        this.dbid = dbid;
        species = org;
    }

    public int getDBID() { return dbid;}    
    public Set<Experiment> getInputs() {return inputs;}
    public Map<String,String> getParams() {return params;}
    public String getParam(String key) {return params.get(key);}
    public Organism getSpecies() {return species;}

}

