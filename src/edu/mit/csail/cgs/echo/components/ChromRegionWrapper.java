/*
 * Created on Feb 19, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.mit.csail.cgs.echo.components;

import java.util.*;

import edu.mit.csail.cgs.ewok.types.ClassType;
import edu.mit.csail.cgs.ewok.types.EchoType;
import edu.mit.csail.cgs.ewok.types.SelfDescribingVerb;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.species.Genome;
import edu.mit.csail.cgs.ewok.nouns.*;

public class ChromRegionWrapper implements SelfDescribingVerb, Generator {
    
    private static EchoType[] paramClasses = { new ClassType(Genome.class) };
    private static String[] paramNames = { "Genome" };
    
    private Map<String,Object> params;
    
    public ChromRegionWrapper() { 
        params = new HashMap<String,Object>();
    }

    public EchoType getOutputClass() {
        return new ClassType(NamedRegion.class);
    }

    public EchoType[] getParameterClasses() { return paramClasses; }
    public String[] getParameterNames() { return paramNames; }
    
    public EchoType[] getInputClasses() { return null; }
    public String[] getInputNames() { return null; }

    public void init(Map<String, Object> params) {
        this.params = new HashMap<String,Object>(params);
    }
    
    public Iterator execute() { 
        Genome g = (Genome)params.get("Genome");
        ChromRegionIterator itr = new ChromRegionIterator(g);
        LinkedList<NamedRegion> regions = new LinkedList<NamedRegion>();
        while(itr.hasNext()) { regions.addLast(itr.next()); }
        
        return regions.iterator();
    }
}
