package edu.mit.csail.cgs.datasets.general;

import edu.mit.csail.cgs.datasets.species.Genome;

public class NamedPoint extends Point implements Named {
    private String name;
    
    public NamedPoint(NamedPoint np) { 
        super(np);
        name = np.name;
    }
    
    public NamedPoint(Point r) { 
        super(r);
        name = r.getLocationString();
        if(r instanceof NamedPoint) { 
            name = ((NamedPoint)r).name;
        }
    }
    
    public NamedPoint(Point r, String n) {
        super(r);
        name = n;
    }
    
    public NamedPoint(Genome g, String c, int loc, String name) {
        super(g,c,loc);
        this.name = name;
    }
        
    public void setName(String n) { name = n; }

    public String getName() {return name;}
    public String toString() {return getName();}
    
    public boolean equals(Object o) { 
        if(!(o instanceof NamedPoint)) { return false; }
        NamedPoint np = (NamedPoint)o;
        if(!name.equals(np.name)) { return false; }
        return super.equals(o);
    }
    
    public int hashCode() { 
        int code = super.hashCode();
        code += name.hashCode(); code *= 37;
        return code; 
    }
}
