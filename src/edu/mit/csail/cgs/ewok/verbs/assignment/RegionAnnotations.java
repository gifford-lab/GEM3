/*
 * Created on Dec 4, 2006
 */
package edu.mit.csail.cgs.ewok.verbs.assignment;

import java.util.*;

import edu.mit.csail.cgs.ewok.nouns.*;
import edu.mit.csail.cgs.ewok.verbs.*;
import edu.mit.csail.cgs.ewok.verbs.binding.*;
import edu.mit.csail.cgs.ewok.verbs.chippet.ChipPetExpander;
import edu.mit.csail.cgs.datasets.binding.*;
import edu.mit.csail.cgs.datasets.chippet.ChipPetDatum;
import edu.mit.csail.cgs.datasets.chippet.ChipPetExpt;
import edu.mit.csail.cgs.datasets.chippet.ChipPetLoader;
import edu.mit.csail.cgs.datasets.general.NamedRegion;
import edu.mit.csail.cgs.datasets.general.Region;
import edu.mit.csail.cgs.datasets.species.Gene;
import edu.mit.csail.cgs.datasets.species.Genome;

/* RegionAnnotations annotations the set of genes in a genome with ChipPet binding events orBindingScan results.
   The set of genes can either be the default from RefGeneGenerator or those from a user-provided
   Gene generator
*/

public class RegionAnnotations extends CachedAnnotations<Region,Region> {
	
    public RegionAnnotations(Genome g) { 
        super();
        
        RefGeneGenerator<NamedRegion> gen = new RefGeneGenerator<NamedRegion>(g);
        Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
        Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
        
        TreeSet<Region> targetList = new TreeSet<Region>();
        while(genes.hasNext()) { 
            Gene gene = genes.next();
            if(!targetList.contains(gene)) { 
                targetList.add(gene);
            }
        }
        
        addItems(targetList);
    }
    
    public RegionAnnotations(Genome g, Expander<NamedRegion,Gene> gen) { 
        super();
        
        Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
        Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
        
        TreeSet<Region> targetList = new TreeSet<Region>();
        while(genes.hasNext()) { 
            Gene gene = genes.next();
            if(!targetList.contains(gene)) { 
                targetList.add(gene);
            }
        }
        
        addItems(targetList);
    }
    
    public RegionAnnotations(Genome g, Expander<NamedRegion,Gene> gen, Mapper<Region,Region> zoneMapper) { 
        super();
        
        Iterator<NamedRegion> chroms = new ChromRegionIterator(g);
        Iterator<Gene> genes = new ExpanderIterator<NamedRegion,Gene>(gen, chroms);
        
        TreeSet<Region> targetList = new TreeSet<Region>();
        while(genes.hasNext()) { 
            Gene gene = genes.next();
            if(!targetList.contains(gene)) { 
                targetList.add(gene);
            }
        }
        
        addItems(targetList);
    }
    
    private class AssignmentExpander<X extends Region> implements Expander<Region,Region> {
    	
    	private AssignmentPredicate pred;
    	private Expander<Region,X> exp;
    	private Mapper<Region,Region> queryMapper;
    	
    	public AssignmentExpander(AssignmentPredicate ap, Expander<Region,X> e) { 
    		pred = ap;
    		exp = e;
    		queryMapper = pred.assignmentZoneMapper();
    	}
    	
    	public Iterator<Region> execute(Region r) { 
    		Region query = queryMapper.execute(r);
    		Iterator<X> events = exp.execute(query);
    		Filter<X,Region> filter = new AssignmentPredicate.Filter(pred, r);
    		return new FilterIterator<X,Region>(filter, events);
    	}
    }
    
    public void addChipPetAnnotations(ChipPetLoader loader, 
    		ChipPetExpt expt, int thresh, AssignmentPredicate ap) {
    	
        ChipPetExpander intvExpander = new ChipPetExpander(loader, expt);
        String key = expt.getName();
        
        addAnnotations(key, new AssignmentExpander(ap, intvExpander));
    }
    
    public void addBindingScanAnnotations(BindingScanLoader loader, 
    		BindingScan scan, AssignmentPredicate ap) {
        
    	BindingScanGenerator scanGen = new BindingScanGenerator(loader, scan);
        String key = scan.toString();
                
        addAnnotations(key, new AssignmentExpander(ap, scanGen));
    }

    public void addBindingScanAnnotations(BindingScanLoader loader, 
    		BindingScan scan, String key, AssignmentPredicate ap) {
        BindingScanGenerator scanGen = new BindingScanGenerator(loader, scan);
        
        addAnnotations(key, new AssignmentExpander(ap, scanGen));
    }
}
